#!/usr/bin/env python
"""Join AF2 scores with peptide_focus into one unfiltered metrics table.

Every design that was predicted gets a row, pass or fail. Thresholds are
applied downstream in funnel.py so that pass rates and filter orthogonality
can be recomputed, and cutoffs swept, without rerunning any GPU work.

Two details about AF2's output that silently corrupt this if missed:

  * It RENUMBERS chains continuously, so the residue-index gap used to mark
    the peptide is gone by output time. The peptide is recovered as the last
    9 residues of chain B instead.
  * Chain A is the binder, chain B the merged HLA + peptide target.

usage: collect_metrics.py <scorefile> <af2_pdb_dir> <out.tsv>
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pmhc_pilot.interface import load_structure, residue_sasa  # noqa: E402

PEPTIDE_LEN = 9


def parse_scorefile(path: Path) -> dict[str, dict]:
    """Parse the Rosetta-style SCORE: table emitted by predict.py."""
    rows: dict[str, dict] = {}
    header: list[str] | None = None
    for line in path.read_text().splitlines():
        if not line.startswith("SCORE:"):
            continue
        fields = line.split()[1:]
        if header is None:
            header = fields
            continue
        if len(fields) != len(header):
            continue
        rec = dict(zip(header, fields))
        tag = rec.get("description")
        if not tag:
            continue
        out: dict[str, float | str] = {"description": tag}
        for key in ("binder_aligned_rmsd", "pae_interaction", "plddt_binder",
                    "plddt_target", "pae_binder", "target_aligned_rmsd"):
            try:
                out[key] = float(rec[key])
            except (KeyError, ValueError):
                out[key] = float("nan")
        rows[tag] = out
    return rows


def peptide_focus(pdb: Path) -> tuple[float, float, float]:
    """Buried peptide area, buried MHC area, and their ratio."""
    model = load_structure(str(pdb))
    chains = {c.id for c in model}
    if not {"A", "B"} <= chains:
        raise ValueError(f"expected chains A and B, found {sorted(chains)}")

    b_res = sorted(
        r.id[1] for c in model if c.id == "B" for r in c if r.id[0] == " "
    )
    if len(b_res) <= PEPTIDE_LEN:
        raise ValueError(f"chain B has only {len(b_res)} residues")
    peptide_ids = set(b_res[-PEPTIDE_LEN:])

    free = residue_sasa(model, ["B"])
    bound = residue_sasa(model, ["A", "B"])

    buried_pep = buried_mhc = 0.0
    for key, area in free.items():
        delta = area - bound.get(key, area)
        if delta <= 0:
            continue
        if key[1] in peptide_ids:
            buried_pep += delta
        else:
            buried_mhc += delta

    total = buried_pep + buried_mhc
    return buried_pep, buried_mhc, (buried_pep / total if total else 0.0)


def backbone_of(tag: str) -> str:
    """pilot_37_dldesign_2_af2pred -> pilot_37, so per-backbone rates work."""
    for marker in ("_dldesign", "_af2pred"):
        if marker in tag:
            tag = tag.split(marker)[0]
    return tag


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    scorefile, pdbdir, out_path = (Path(a) for a in sys.argv[1:])

    scores = parse_scorefile(scorefile)
    print(f"  {len(scores)} rows in scorefile")

    columns = ["backbone_id", "seq_id", "description", "scrmsd",
               "plddt_binder", "pae_interaction", "peptide_focus",
               "buried_peptide", "buried_mhc", "plddt_target",
               "target_aligned_rmsd"]

    written = failed = 0
    with out_path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for tag, rec in sorted(scores.items()):
            pdb = pdbdir / f"{tag}.pdb"
            try:
                bp, bm, focus = peptide_focus(pdb)
            except Exception as exc:
                # Record the row with a missing focus rather than dropping it:
                # a silently absent design inflates every pass rate.
                bp = bm = focus = float("nan")
                failed += 1
                if failed <= 3:
                    print(f"  focus failed for {tag}: {exc}")
            fh.write("\t".join(str(v) for v in [
                backbone_of(tag), tag.split("_dldesign_")[-1].split("_")[0],
                tag, rec["binder_aligned_rmsd"], rec["plddt_binder"],
                rec["pae_interaction"], f"{focus:.4f}", f"{bp:.1f}",
                f"{bm:.1f}", rec["plddt_target"], rec["target_aligned_rmsd"],
            ]) + "\n")
            written += 1

    print(f"  wrote {written} rows ({failed} without focus) -> {out_path}")


if __name__ == "__main__":
    main()
