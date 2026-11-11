#!/usr/bin/env python
"""Convert RFdiffusion output into the two-chain form dl_binder_design wants.

dl_binder_design requires exactly two chains, binder first:

    predict.py:  if len(splits) > 2: raise Exception(...)
    dl_interface_design.py:  fixed_res = {chains[0]: ..., chains[1]: []}

RFdiffusion's two modes emit DIFFERENT layouts, and both need work:

  free binder design (run_rfdiffusion.sh)
      A = HLA platform (1-180), B = binder (1-130), C = peptide (1-9)
      Three chains; must be merged down to two.

  scaffold-guided / fold-conditioned (run_rfdiffusion_scaffolded.sh)
      A = binder (1-N), B = HLA + peptide (N+1 .. N+189), CONTIGUOUS
      Already two chains, but the HLA and peptide are fused into one
      unbroken numbering -- so the peptide reads as covalently bonded to
      the HLA C-terminus, which it is not.

Both are normalised to the same output: chain A = binder numbered from 1,
chain B = HLA numbered from 1 with the peptide renumbered after a large gap
(PEPTIDE_OFFSET). A residue-index discontinuity is how AF2 represents a
chain break, so the model does not try to bridge them.

The target is held fixed by the initial-guess protocol in any case, but a
spurious bond would still distort the interface we are trying to measure.

Emitting one format from both modes is what lets the rest of the pipeline --
collect_metrics.py in particular, which recovers the peptide as the last
nine residues of chain B -- stay identical across arms. A comparison between
arms is only as trustworthy as the sameness of everything downstream of the
variable being compared.
"""

import os
import sys
from pathlib import Path

# Residue-index gap inserted between the MHC's C-terminus and the peptide, to
# read as a chain break rather than a linker.
#
# This is a GAP, not an absolute offset. It used to be `PEPTIDE_OFFSET = 200`,
# with the peptide renumbered from 201 -- but the MHC is 180 residues, so the
# resulting gap was 21, which sits INSIDE the +/-32 window AF2's relative
# position encoding represents. A 21-residue separation reads as a plausible
# loop, not as two molecules. Computing the gap from the MHC's actual length
# also keeps this correct for a target of a different size, which the transfer
# arm needs.
PEPTIDE_GAP = 64  # > 32, so it saturates the clip whatever the MHC length

# The 9-mer in the groove. Used to verify the MHC/peptide split in
# scaffold-guided output, where the boundary is positional rather than
# marked by a chain id -- a silent off-by-one there would misattribute
# buried surface between the numerator and denominator of peptide_focus.
# Overridable so a second target (KF9 on HLA-C*12:02) can reuse this
# pipeline. Default unchanged, so every NY-ESO-1 arm behaves exactly as
# before and the assertion still fires on a wrong peptide boundary.
EXPECTED_PEPTIDE = os.environ.get("PMHC_PEPTIDE", "SLLMWITQV")

# 9MIN trimmed to the HLA a1/a2 platform by scripts/prep_target.py. Used to
# tell the MHC chain from the binder when chain letters are unreliable; an
# absolute length, because a RELATIVE test ("the binder is the smaller one")
# is tautological once the chains are sorted by size and can never fail.
MHC_LEN = 180

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _by_chain(src: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in src.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out.setdefault(line[21], []).append(line)
    return out


def _residues(lines: list[str]) -> list[tuple[str, list[str]]]:
    """Group atom lines into residues, preserving file order."""
    groups: list[tuple[str, list[str]]] = []
    index: dict[str, list[str]] = {}
    for ln in lines:
        key = ln[22:27]  # resseq + icode
        if key not in index:
            index[key] = []
            groups.append((key, index[key]))
        index[key].append(ln)
    return groups


def _sequence(residues) -> str:
    return "".join(THREE_TO_ONE.get(atoms[0][17:20].strip(), "X")
                   for _, atoms in residues)


def _emit(out: list[str], counts: dict[str, int], residues, chain: str,
          start: int) -> None:
    """Append residues under `chain`, renumbered consecutively from `start`."""
    for offset, (_, atoms) in enumerate(residues):
        num = start + offset
        for ln in atoms:
            out.append(f"{ln[:21]}{chain}{num:>4} {ln[27:]}")
            if ln[12:16].strip() == "CA":
                counts[chain] = counts.get(chain, 0) + 1


def _write(dst: Path, binder, mhc, peptide) -> dict[str, int]:
    out: list[str] = []
    counts: dict[str, int] = {}
    # Binder first -- the pipeline predicts chains[0] and fixes the rest.
    _emit(out, counts, binder, "A", 1)
    _emit(out, counts, mhc, "B", 1)
    _emit(out, counts, peptide, "B", len(mhc) + PEPTIDE_GAP)
    out += ["TER", "END"]
    dst.write_text("\n".join(out) + "\n")
    return counts


def convert_three_chain(src: Path, dst: Path,
                        peptide_len=len(EXPECTED_PEPTIDE),
                        mhc_len=MHC_LEN) -> dict[str, int]:
    """Three-chain output -> two chains, identifying roles by size.

    The chain LETTERS are not stable across RFdiffusion modes:

        free binder design   A = MHC (180), B = binder, C = peptide (9)
        partial diffusion    A = binder, B = MHC (180), C = peptide (9)

    Dispatching on the letter would silently treat a 128-residue binder as
    the MHC in the second case, which does not error -- it just computes
    peptide_focus against the wrong chains. So roles are assigned by size
    instead: the peptide is the chain of exactly `peptide_len`, the MHC is
    the largest of the rest, and the binder is what remains.
    """
    chains = _by_chain(src)
    if len(chains) != 3:
        raise ValueError(f"expected 3 chains, found {sorted(chains)}")

    grouped = {cid: _residues(lines) for cid, lines in chains.items()}
    sizes = {c: len(r) for c, r in grouped.items()}

    def exactly_one(length: int, role: str) -> str:
        hits = [c for c, n in sizes.items() if n == length]
        if len(hits) != 1:
            raise ValueError(
                f"need exactly one {length}-residue chain to be the {role}, "
                f"found {len(hits)}; chain sizes are {sizes}")
        return hits[0]

    pep = exactly_one(peptide_len, "peptide")
    mhc = exactly_one(mhc_len, "MHC")
    rest = [c for c in grouped if c not in (pep, mhc)]
    if len(rest) != 1:
        raise ValueError(f"cannot identify the binder chain; sizes are {sizes}")
    return _write(dst, grouped[rest[0]], grouped[mhc], grouped[pep])


# Retained for callers that know the layout; dispatch prefers the size rule.
def convert_free(src: Path, dst: Path, binder_src="B", mhc_src="A",
                 peptide_src="C") -> dict[str, int]:
    """Three-chain free-design output -> two chains, by explicit chain ids."""
    chains = _by_chain(src)
    missing = {binder_src, mhc_src, peptide_src} - set(chains)
    if missing:
        raise ValueError(f"missing chains {sorted(missing)}")
    return _write(dst,
                  _residues(chains[binder_src]),
                  _residues(chains[mhc_src]),
                  _residues(chains[peptide_src]))


def convert_scaffolded(src: Path, dst: Path, binder_src="A", target_src="B",
                       peptide_len=len(EXPECTED_PEPTIDE)) -> dict[str, int]:
    """Scaffold-guided output -> two chains, with the chain break restored.

    The peptide is the tail of the fused target chain. Its identity is
    checked rather than assumed: the split point carries no chain id, so
    nothing but the sequence would catch an off-by-one.
    """
    chains = _by_chain(src)
    missing = {binder_src, target_src} - set(chains)
    if missing:
        raise ValueError(f"missing chains {sorted(missing)}")

    target = _residues(chains[target_src])
    if len(target) <= peptide_len:
        raise ValueError(f"target chain has only {len(target)} residues")
    mhc, peptide = target[:-peptide_len], target[-peptide_len:]

    seq = _sequence(peptide)
    if seq != EXPECTED_PEPTIDE:
        raise ValueError(
            f"tail of chain {target_src} is {seq!r}, expected "
            f"{EXPECTED_PEPTIDE!r} -- the peptide boundary is wrong, so "
            f"peptide_focus would be computed against the wrong residues"
        )
    return _write(dst, _residues(chains[binder_src]), mhc, peptide)


def convert(src: Path, dst: Path) -> dict[str, int]:
    """Dispatch on chain layout so one pipeline handles every arm."""
    ids = set(_by_chain(src))
    if len(ids) == 3:
        return convert_three_chain(src, dst)
    if len(ids) == 2:
        return convert_scaffolded(src, dst)
    raise ValueError(f"unrecognised chain layout {sorted(ids)}")


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit("usage: to_two_chain.py <in_dir_or_pdb> <out_dir>")

    src_arg = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    sources = sorted(src_arg.glob("*.pdb")) if src_arg.is_dir() else [src_arg]
    if not sources:
        sys.exit(f"no PDBs found in {src_arg}")

    converted = failed = 0
    for src in sources:
        try:
            counts = convert(src, out_dir / src.name)
            converted += 1
            if converted <= 3:
                print(f"  {src.name}: A={counts.get('A')} B={counts.get('B')}")
        except Exception as exc:
            failed += 1
            print(f"  {src.name}: FAILED {exc}")

    print(f"\nconverted {converted}, failed {failed} -> {out_dir}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
