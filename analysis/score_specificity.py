#!/usr/bin/env python
"""Fifth filter: does the designed interface PREFER the on-target peptide?

    usage:  score_specificity.py <metrics.tsv> <af2_pdb_dir> <out.tsv>
                                 [--peptide SLLMWITQV] [--samples 10]
                                 [--all] [--decoys FILE]

WHAT THIS MEASURES AND WHY THE FOUR EXISTING GATES CANNOT

scRMSD, pLDDT, interface pAE and peptide_focus all ask "is this a confident
complex?". None asks "is it confident about THIS peptide rather than any
peptide in the groove". A binder that clamps the HLA helices and ignores the
9-mer scores well on all four. That is the single largest gap between this
pilot and the published pipeline, and it is the cheapest one to close,
because it runs on designs that are already scored and needs no new GPU work.

THE METHOD, AND WHY IT IS ONE PASS RATHER THAN NINE

Liu et al.: "ProteinMPNN was used to score the likelihood of each peptide
residue being predicted when in complex with a binder. Binders were screened
by scoring each peptide residue separately and filtering directly on the log
probability of the on-target residue, or on the difference in log
probabilities to when that residue is mutated to an alanine."

ProteinMPNN's --conditional_probs_only returns exactly p(s_i | s_{j!=i}, X)
for every position in a single forward pass -- each residue already scored
with the other eight held fixed. It returns the full 21-way distribution, so
log P(native) and log P(Ala) at the same position come from the same numbers
and the alanine-difference filter costs nothing extra. Nine separate runs
would compute the same thing nine times more slowly.

MPNN is BACKBONE-ONLY. It never sees a side chain, of the peptide, the
binder or the MHC. So "mutating" a peptide residue means relabelling it; no
repacking, no threading, no structural work at all. That is what makes an
off-target panel affordable (--decoys).

WHY holo - apo RATHER THAN holo ALONE

Raw log P(peptide | complex) is dominated by the groove. HLA-A*02:01 has
strong anchor preferences -- Leu at P2, Val at P9 -- and a design that reads
nothing still inherits them, so the absolute number mostly measures the MHC.
Subtracting the same quantity computed with the binder chain DELETED removes
that baseline and leaves the binder's own contribution:

    delta_i = log P(s_i | binder present) - log P(s_i | binder absent)

Positive means the binder makes its target residue more predictable than the
groove alone does, which is the thing we actually want to filter on. Measured
on one NY-ESO-1 design: raw holo was -25.5 against a -27.0 uniform baseline,
nearly useless; holo-apo was +2.74 on-target against -1.90 for a KRAS decoy.

SCOPE

By default only designs passing all four AF2 gates are scored, because a
specificity number for a design that is not a credible complex means nothing.
--all overrides this.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
ALA = MPNN_ALPHABET.index("A")

# Same four gates as pmhc_pilot.funnel, AF2-shifted. Kept explicit so this
# script can run against a metrics file without importing the package.
GATES = (
    ("scrmsd",          lambda v: v < 2.0),
    ("plddt_binder",    lambda v: v > 92.0),
    ("pae_interaction", lambda v: v < 7.0),
    ("peptide_focus",   lambda v: v >= 0.202),
)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}


def read_pdb_chains(path):
    """chain -> ordered list of (resnum, [lines])."""
    chains, order = {}, {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[16] not in (" ", "A"):
            continue
        ch, num = line[21], int(line[22:26])
        per = chains.setdefault(ch, {})
        if num not in per:
            per[num] = []
            order.setdefault(ch, []).append(num)
        per[num].append(line)
    return {c: [(n, chains[c][n]) for n in order[c]] for c in order}


def write_variant(src_chains, dst, keep, peptide_seq=None, pep_len=9):
    """Write selected chains; optionally relabel the trailing peptide.

    Relabelling only rewrites the residue NAME columns. Backbone coordinates
    are untouched, which is all ProteinMPNN reads -- so a decoy peptide costs
    no structural modelling.
    """
    out = []
    for ch in keep:
        residues = src_chains[ch]
        for idx, (_, lines) in enumerate(residues):
            new_name = None
            if peptide_seq is not None and ch == keep[-1]:
                offset = idx - (len(residues) - pep_len)
                if offset >= 0:
                    new_name = ONE_TO_THREE[peptide_seq[offset]]
            for ln in lines:
                out.append(ln[:17] + (new_name or ln[17:20]) + ln[20:])
        out.append("TER\n")
    out.append("END\n")
    Path(dst).write_text("".join(out))


def run_mpnn(mpnn_root, pdb_dir, out_dir, samples):
    """One batched conditional-probs run over a directory of PDBs."""
    jsonl = Path(out_dir) / "in.jsonl"
    subprocess.run(
        [sys.executable, str(Path(mpnn_root) / "helper_scripts"
                             / "parse_multiple_chains.py"),
         "--input_path", str(pdb_dir), "--output_path", str(jsonl)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(
        [sys.executable, str(Path(mpnn_root) / "protein_mpnn_run.py"),
         "--jsonl_path", str(jsonl), "--out_folder", str(out_dir),
         "--conditional_probs_only", "1",
         "--num_seq_per_target", str(samples), "--batch_size", "1"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    probs = {}
    for npz in (Path(out_dir) / "conditional_probs_only").glob("*.npz"):
        d = np.load(npz)
        # log_p is [num_decoding_orders, L, 21]; average in probability space,
        # not log space -- these are samples of the same distribution under
        # different decoding orders, and averaging logs would understate it.
        probs[npz.stem] = (np.log(np.exp(d["log_p"]).mean(axis=0) + 1e-12),
                           d["S"])
    return probs


def passes_gates(row):
    try:
        return all(test(float(row[col])) for col, test in GATES)
    except (KeyError, ValueError, TypeError):
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metrics")
    ap.add_argument("pdb_dir")
    ap.add_argument("out")
    ap.add_argument("--peptide", default="SLLMWITQV")
    ap.add_argument("--samples", type=int, default=10,
                    help="decoding orders to average over (variance control)")
    ap.add_argument("--all", action="store_true",
                    help="score every design, not just four-gate survivors")
    ap.add_argument("--decoys",
                    help="file of off-target peptides, one per line, same length")
    ap.add_argument("--mpnn-root",
                    default=str(Path.home() / "pmhc/dl_binder_design/mpnn_fr/ProteinMPNN"))
    args = ap.parse_args()

    pep = args.peptide.upper()
    decoys = []
    if args.decoys:
        decoys = [l.strip().upper() for l in open(args.decoys)
                  if l.strip() and not l.startswith("#")]
        bad = [d for d in decoys if len(d) != len(pep)]
        if bad:
            sys.exit(f"decoys must be {len(pep)}-mers; offenders: {bad[:3]}")

    import csv
    rows = list(csv.DictReader(open(args.metrics), delimiter="\t"))
    if not args.all:
        rows = [r for r in rows if passes_gates(r)]
    if not rows:
        sys.exit("no designs to score (no four-gate survivors; use --all)")

    pdb_dir = Path(args.pdb_dir)
    selected = []
    for r in rows:
        hit = pdb_dir / f"{r['description']}.pdb"
        if hit.is_file():
            selected.append((r["description"], hit))
    missing = len(rows) - len(selected)
    print(f"{len(selected)} designs to score"
          + (f" ({missing} had no PDB in {pdb_dir})" if missing else ""))
    if not selected:
        sys.exit(f"none of the designs were found in {pdb_dir}")

    work = Path(tempfile.mkdtemp(prefix="spec_"))
    try:
        # Stage holo (binder + pMHC) and apo (pMHC only), plus one staging
        # directory per decoy peptide. Every variant of every design goes
        # through MPNN in a single batched run per condition.
        conditions = {"holo": None, "apo": None}
        for d in decoys:
            conditions[f"decoy:{d}"] = d

        staged, lengths = {}, {}
        for cond in conditions:
            (work / cond.replace(":", "_")).mkdir(parents=True, exist_ok=True)
        for name, path in selected:
            chains = read_pdb_chains(path)
            if "A" not in chains or "B" not in chains:
                continue
            lengths[name] = (len(chains["A"]), len(chains["B"]))
            for cond, decoy in conditions.items():
                tag = cond.replace(":", "_")
                keep = ["B"] if cond == "apo" else ["A", "B"]
                write_variant(chains, work / tag / f"{name}.pdb", keep,
                              peptide_seq=decoy, pep_len=len(pep))
            staged[name] = True

        results = {}
        for cond in conditions:
            tag = cond.replace(":", "_")
            print(f"  MPNN pass: {cond} ({len(staged)} structures)...", flush=True)
            results[cond] = run_mpnn(args.mpnn_root, work / tag,
                                     work / f"out_{tag}", args.samples)

        with open(args.out, "w") as fh:
            fh.write("design\tcondition\tposition\tresidue\tlogp_holo\t"
                     "logp_apo\tdelta\tlogp_ala\tdelta_ala\n")
            summary = {}
            for name in staged:
                if name not in results["holo"] or name not in results["apo"]:
                    continue
                lp_h, S_h = results["holo"][name]
                lp_a, _ = results["apo"][name]
                n = len(pep)
                for cond in conditions:
                    lp_c, S_c = results[cond].get(name, (None, None))
                    if lp_c is None:
                        continue
                    seq = (pep if cond in ("holo", "apo")
                           else cond.split(":", 1)[1])
                    tot = 0.0
                    for i in range(n):
                        h_idx, a_idx = -n + i, -n + i      # trailing peptide
                        aa = MPNN_ALPHABET.index(seq[i])
                        # sanity: the staged sequence must be what MPNN parsed
                        if S_c[h_idx] != aa:
                            print(f"    WARNING {name} {cond} pos {i+1}: "
                                  f"expected {seq[i]}, MPNN read "
                                  f"{MPNN_ALPHABET[S_c[h_idx]]}")
                        lh = float(lp_c[h_idx, aa])
                        la = float(lp_a[a_idx, aa])
                        lala = float(lp_c[h_idx, ALA])
                        tot += lh - la
                        fh.write(f"{name}\t{cond}\t{i+1}\t{seq[i]}\t"
                                 f"{lh:.4f}\t{la:.4f}\t{lh-la:.4f}\t"
                                 f"{lala:.4f}\t{lh-lala:.4f}\n")
                    summary.setdefault(name, {})[cond] = tot

        print(f"\nwrote {args.out}")
        on = [v["holo"] for v in summary.values() if "holo" in v]
        if on:
            print(f"  summed holo-apo, on-target: median {np.median(on):+.2f}  "
                  f"range {min(on):+.2f} to {max(on):+.2f}  (n={len(on)})")
        for d in decoys:
            k = f"decoy:{d}"
            vals = [v[k] for v in summary.values() if k in v]
            if vals:
                better = sum(1 for nm, v in summary.items()
                             if k in v and "holo" in v and v["holo"] > v[k])
                print(f"  decoy {d}: median {np.median(vals):+.2f}   "
                      f"on-target wins in {better}/{len(vals)}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
