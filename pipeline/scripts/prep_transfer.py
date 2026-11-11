#!/usr/bin/env python
"""Build the partial-diffusion transfer input: KH46 placed on a NEW pMHC.

The transfer arm asks whether a validated pMHC-binding ARCHITECTURE carries
over to a different peptide in the same groove. Partial diffusion starts from
coordinates rather than from a topology abstraction, so the binder has to
already be sitting on the new target before diffusion begins. That is what
this script constructs.

THE METHOD

Both structures are HLA-A*02:01 -- verified, not assumed: the alpha1/alpha2
platforms are 180/180 IDENTICAL in sequence, so this is the same molecule in
two crystal frames (0.403 A CA RMSD). Superimposing the new platform onto
9MIN's and carrying the peptide along with it places the new peptide in the
groove exactly where the old one sat, under an unmoved KH46.

The peptide is NOT threaded. Its sequence is not modelled onto 9MIN's peptide
backbone; the whole 2GUO peptide is moved as a rigid body, so it keeps the
conformation it was crystallised in. This matters because the two epitopes do
not adopt the same backbone: SLLMWITQV presents a prominent Trp5, MART-1 is
comparatively flat. Threading would have forced MART-1 into a conformation it
does not adopt and the cavity reported below would have been an artefact of
the modelling rather than a real difference between the two antigens.

TWO MODES, BECAUSE "SAME TARGET" IS AMBIGUOUS

    graft  (default)  9MIN's own MHC + the moved 2GUO peptide
    whole             2GUO's MHC + 2GUO's peptide, both moved

`graft` makes the peptide the only changed thing: the MHC surface is bit
identical to the one KH46 was designed against, and identical to the control
arm's target, so any change in MHC contacts is the peptide reshaping what the
groove presents rather than crystal-form noise. It is also cleaner -- KH46
against 9MIN's MHC is its own crystal structure, so zero clashes there by
construction.

`whole` keeps the genuine MART-1-bound groove conformation, which `graft`
gives up: 9MIN's groove is adapted to SLLMWITQV. That is a real objection,
but a second-order one at 0.403 A, and it costs the clean single variable. The binder keeps its
original coordinates; the target moves underneath it.

WHICH STRUCTURE, AND THE TRAP

"MART-1" names at least three different peptides in the PDB, and two of the
obvious candidates are the wrong length or the wrong sequence:

    1JF1   ELAGIGILTV   10mer, A27L anchor-modified   WRONG LENGTH
    1JHT   ALGIGILTV     9mer, anchor-modified        not the natural epitope
    2GUO   AAGIGILTV     9mer, natural MART-1 27-35   <-- this one
    3QEU   (no MHC)      free DMF5 TCR                unusable

A decamer would have been silently wrong rather than loudly wrong: every
length assumption downstream is 9 (to_two_chain's peptide_len, the chain-break
offset, collect_metrics recovering the peptide as the last nine residues of
chain B). The nonamer means the transfer arm needs no code changes at all,
which is the whole point -- the arm has to differ from the others in the
target and nothing else.

2GUO also carries two copies in the asymmetric unit; the second (D/E/F) is
dropped, as is beta-2-microglobulin and the alpha3 domain, exactly as
prep_target.py does for 9MIN.

WHAT THIS DOES NOT DO

It does not relax or repack the interface. The binder is placed and left
alone, so whatever KH46 does against the new peptide -- clash, cavity, or
fit -- is reported honestly and handed to diffusion to resolve. Measuring
that gap before spending GPU time is the point of the clash report below.
"""

import sys
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "targets"

# The binder and the 9MIN reference platform, as prep_target.py wrote them.
REFERENCE = OUT / "9min_control_kh46.pdb"
REF_MHC, REF_PEP, REF_BINDER = "A", "C", "D"

# alpha1 + alpha2; alpha3 begins ~183. Must match prep_target.PLATFORM_LAST,
# because to_two_chain identifies the MHC chain by this exact length.
PLATFORM_LAST = 180

# Defined here rather than imported from to_two_chain: this script is run with
# python -I, which drops the script's own directory from sys.path.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# A contact closer than this between two heavy atoms is a steric clash.
CLASH = 2.5
# Heavy atoms within this of each other are in contact; used to show that the
# transferred binder still engages the groove at all.
CONTACT = 4.5


def fetch(stem: str) -> Path:
    """Ensure data/raw/<stem>.cif exists, downloading it if not.

    data/raw is gitignored, so a fresh clone has no structures. Fetching here
    rather than in a separate setup step keeps this script runnable on its own.
    """
    path = RAW / f"{stem}.cif"
    if path.exists():
        return path
    RAW.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{stem.upper()}.cif"
    print(f"fetching {url}")
    urllib.request.urlretrieve(url, path)
    return path


def heavy(residues):
    """Backbone + side-chain heavy atoms, dropping hydrogens and altlocs."""
    out = []
    for r in residues:
        for a in r:
            if a.element != "H" and a.get_altloc() in (" ", "A"):
                out.append(a)
    return out


def ordered(chain, last=None):
    rs = [r for r in chain if r.id[0] == " "]
    return [r for r in rs if last is None or r.id[1] <= last]


def emit(handle, residues, chain_id, start, serial):
    """Write residues as chain_id, renumbered consecutively from start."""
    for offset, r in enumerate(residues):
        for a in r:
            if a.element == "H" or a.get_altloc() not in (" ", "A"):
                continue
            x, y, z = a.coord
            name = a.get_name()
            # PDB atom names are left-justified from column 14 unless the
            # element symbol is two characters.
            nm = f" {name:<3}" if len(name) < 4 else name
            handle.write(
                f"ATOM  {serial:>5} {nm}{' '}{r.get_resname():>3} "
                f"{chain_id}{start + offset:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                f"{a.get_bfactor():>6.2f}{0.0:>6.2f}\n"
            )
            serial += 1
    handle.write("TER\n")
    return serial


def report_interface(binder, target, label):
    """Min distance, clash count and contact count, binder vs target."""
    b = np.array([a.coord for a in heavy(binder)])
    t = np.array([a.coord for a in heavy(target)])
    d = np.linalg.norm(b[:, None, :] - t[None, :, :], axis=-1)
    print(f"  {label:<22} min {d.min():>6.2f} A   "
          f"clashes <{CLASH} {int((d < CLASH).sum()):>4}   "
          f"contacts <{CONTACT} {int((d < CONTACT).sum()):>5}")
    return d.min(), int((d < CLASH).sum()), int((d < CONTACT).sum())


def main() -> None:
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = {a for a in sys.argv[1:] if a.startswith("--")}
    mode = "whole" if "--whole" in flags else "graft"
    if flags - {"--whole", "--graft"}:
        sys.exit(f"unknown flag(s) {sorted(flags - {'--whole', '--graft'})}; "
                 f"usage: prep_transfer.py [pdb_id] [mhc_chain] [pep_chain] "
                 f"[--graft|--whole]")

    stem = (argv[0] if argv else "2guo").lower()
    stem = stem[:-4] if stem.endswith(".cif") else stem
    cif = fetch(stem)
    mhc_src, pep_src = (argv[1], argv[2]) if len(argv) > 2 else ("A", "C")

    ref = PDBParser(QUIET=True).get_structure("ref", REFERENCE)[0]
    new = MMCIFParser(QUIET=True).get_structure("new", cif)[0]

    ref_mhc = ordered(ref[REF_MHC], PLATFORM_LAST)
    ref_pep = ordered(ref[REF_PEP])
    binder = ordered(ref[REF_BINDER])
    new_mhc = ordered(new[mhc_src], PLATFORM_LAST)
    new_pep = ordered(new[pep_src])

    print(f"reference  {REFERENCE.name}")
    print(f"  MHC {len(ref_mhc)}   peptide {len(ref_pep)}   "
          f"binder {len(binder)}")
    print(f"new target {cif.name}  chains {mhc_src}/{pep_src}")
    print(f"  MHC {len(new_mhc)}   peptide {len(new_pep)}")

    if len(new_mhc) != len(ref_mhc):
        sys.exit(f"platform length mismatch: {len(new_mhc)} vs "
                 f"{len(ref_mhc)}; to_two_chain identifies the MHC by an "
                 f"exact length of {PLATFORM_LAST}")
    if len(new_pep) != len(ref_pep):
        sys.exit(f"peptide length mismatch: {len(new_pep)} vs "
                 f"{len(ref_pep)}; every length assumption downstream is "
                 f"{len(ref_pep)} and a different one would be silently wrong")

    # Superimpose on CA of residues present in BOTH platforms, matched by
    # residue number rather than position, so a missing residue cannot shift
    # the whole alignment by one.
    ref_by_num = {r.id[1]: r for r in ref_mhc}
    new_by_num = {r.id[1]: r for r in new_mhc}
    shared = sorted(set(ref_by_num) & set(new_by_num))

    # SAME ALLELE? Verified by sequence, not assumed from the literature. If
    # the platforms differed, the superposition would be dropping a peptide
    # into a groove whose pockets are not the ones it binds, and the anchor
    # residues would be in the wrong place.
    mism = [(n, THREE_TO_ONE.get(ref_by_num[n].get_resname(), "X"),
             THREE_TO_ONE.get(new_by_num[n].get_resname(), "X"))
            for n in shared
            if ref_by_num[n].get_resname() != new_by_num[n].get_resname()]
    ident = 100.0 * (len(shared) - len(mism)) / len(shared)
    print(f"\nallele check: {len(shared) - len(mism)}/{len(shared)} platform "
          f"residues identical ({ident:.1f}%)")
    if mism:
        print(f"  substitutions: "
              + ", ".join(f"{a}{n}{b}" for n, a, b in mism[:12]))
    if ident < 95.0:
        sys.exit(f"platforms are only {ident:.1f}% identical -- these are "
                 f"different alleles, and the groove pockets differ. A "
                 f"superposition would misplace the peptide's anchors.")

    pairs = [(ref_by_num[n]["CA"], new_by_num[n]["CA"]) for n in shared
             if "CA" in ref_by_num[n] and "CA" in new_by_num[n]]
    print(f"\nsuperposition on {len(pairs)} platform CA atoms")

    sup = Superimposer()
    sup.set_atoms([p[0] for p in pairs], [p[1] for p in pairs])
    print(f"  platform RMSD {sup.rms:.3f} A")
    if sup.rms > 1.5:
        sys.exit(f"platform RMSD {sup.rms:.2f} A is too high for two copies "
                 f"of the same allele; check the chain assignment")

    # Move MHC and peptide TOGETHER, as one rigid body, so the peptide keeps
    # the groove relationship it was crystallised in rather than being
    # threaded onto the reference peptide's backbone.
    sup.apply([a for r in new_mhc + new_pep for a in r])

    # Which MHC the binder actually sees. See the two modes in the docstring.
    out_mhc = ref_mhc if mode == "graft" else new_mhc
    print(f"\nmode: {mode}  -- MHC supplied to diffusion is "
          f"{'9MIN (reference)' if mode == 'graft' else cif.stem.upper()}")

    if mode == "graft":
        # The grafted peptide comes from a frame 0.4 A away, so check it does
        # not collide with the reference groove it is being dropped into.
        print("  grafted peptide vs the reference groove:")
        report_interface(new_pep, ref_mhc, "new peptide / 9MIN MHC")

    print("\ninterface, binder vs target")
    print("  --- original complex, for scale ---")
    report_interface(binder, ref_pep, "KH46 / old peptide")
    report_interface(binder, ref_mhc, "KH46 / old MHC")
    print("  --- after transfer ---")
    pep_min, pep_clash, pep_contact = report_interface(
        binder, new_pep, "KH46 / new peptide")
    report_interface(binder, out_mhc, "KH46 / new MHC")

    dst = OUT / f"transfer_{cif.stem}_kh46_{mode}.pdb"
    with dst.open("w") as fh:
        # A = binder, B = MHC, C = peptide: the layout run_partial_transfer.sh
        # reads chain lengths from, and to_two_chain resolves by size.
        serial = emit(fh, binder, "A", 1, 1)
        serial = emit(fh, out_mhc, "B", 1, serial)
        emit(fh, new_pep, "C", 1, serial)
        fh.write("END\n")

    pep_seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in new_pep)
    print(f"\nwrote {dst.name}")
    print(f"  A={len(binder)} (KH46)  B={len(out_mhc)} (MHC)  "
          f"C={len(new_pep)} ({pep_seq})")

    if pep_clash:
        print(f"\n  NOTE: {pep_clash} clashing atom pairs against the new "
              f"peptide.\n  Partial diffusion must resolve these; a high "
              f"count argues for a larger partial_T.")
    elif pep_contact < 50:
        print(f"\n  NOTE: only {pep_contact} contacts to the new peptide "
              f"(no clashes).\n  The transferred binder under-engages it -- "
              f"a cavity, not a collision.")


if __name__ == "__main__":
    main()
