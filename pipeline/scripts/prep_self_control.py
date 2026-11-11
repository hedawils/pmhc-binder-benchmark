#!/usr/bin/env python3
"""Build the KH46-onto-its-own-target input for the partial-diffusion control.

WHY THIS ARM EXISTS

Arm 1 conditions on KH46's secondary structure and block adjacency. Arms 2 and
3 instead partially diffuse from a docked pose. So arm 1 differs from them in
BOTH method and starting material, and cannot isolate either:

    arm 2  KH46 pose    -> MART-1       partial diffusion
    arm 3  152 scaffolds -> NY-ESO-1    partial diffusion
    arm 1  KH46 ss/adj   -> NY-ESO-1    fold conditioning     <-- not matched

This arm supplies the missing cell. Partially diffusing KH46 onto its OWN
target is method-matched to both, which turns two confounded comparisons into
clean ones:

    vs arm 2   same source, same method, same partial_T sweep; only the
               PEPTIDE changes -> isolates the target, and removes the
               "arm 2 starts from a validated binder" confound
    vs arm 3   same target, same method; only the SOURCE changes ->
               isolates source quality, validated-for-this-peptide against
               validated-for-other-peptides

WHAT IT IS AND IS NOT

It is a positive control: a ceiling and a sanity check. At low partial_T it
will largely regenerate a binder that already works, so a high pass rate here
is the expected result and not evidence of anything. Its value is as a
denominator for the other two arms, not as a demonstration.

WHAT THIS SCRIPT DOES

Nothing but relabel. prep_target.py wrote the reference complex as
A=MHC(180) C=peptide(9) D=binder(128); run_partial_transfer.sh requires
A=binder B=MHC C=peptide. Chains are rewritten into that order and residues
renumbered from 1, with columns 31 onward -- x, y, z, occupancy, B, element --
copied byte for byte. No superposition, no threading, no geometry of any kind.
The coordinate multiset is asserted identical before the file is written, so a
regression here fails loudly rather than silently moving the pose.

    usage:  prep_self_control.py <9min_control_kh46.pdb> <out.pdb>
"""

import sys
from collections import OrderedDict
from pathlib import Path

# prep_target.py's layout -> run_partial_transfer.sh's layout.
#   source chain, destination chain, expected CA count
PLAN = (("D", "A", 128), ("A", "B", 180), ("C", "C", 9))

EXPECTED_PEPTIDE = "SLLMWITQV"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def read_chains(path):
    """chain id -> OrderedDict of (resseq, icode) -> list of ATOM lines."""
    chains = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[16] not in (" ", "A"):      # altloc: keep the first only
                continue
            chains.setdefault(line[21], OrderedDict()) \
                  .setdefault((line[22:26], line[26]), []).append(line)
    return chains


def sequence(residues):
    return "".join(
        THREE_TO_ONE.get(lines[0][17:20].strip(), "X")
        for lines in residues.values()
    )


def coord_set(lines):
    """Multiset of the coordinate columns, as written."""
    return sorted(ln[30:54] for ln in lines)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__.strip().splitlines()[-1].strip())
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    if not src.is_file():
        sys.exit(f"no such input: {src}")

    chains = read_chains(src)
    print(f"{src.name}: chains " +
          "  ".join(f"{c}={len(r)}" for c, r in chains.items()))

    for source, _, expect in PLAN:
        if source not in chains:
            sys.exit(f"input has no chain {source}; expected "
                     + " ".join(f"{s}({n})" for s, _, n in PLAN))
        got = len(chains[source])
        if got != expect:
            sys.exit(f"chain {source} has {got} residues, expected {expect}. "
                     "This script is pinned to prep_target.py's layout; if "
                     "that changed, fix PLAN rather than loosening the check.")

    pep = sequence(chains["C"])
    if pep != EXPECTED_PEPTIDE:
        sys.exit(f"chain C reads {pep}, expected {EXPECTED_PEPTIDE}. This arm "
                 "is the NY-ESO-1 self-control; a different peptide means the "
                 "wrong input file.")
    print(f"  peptide {pep}  (NY-ESO-1 157-165, as expected)")

    out, serial = [], 1
    for source, destination, _ in PLAN:
        residues = chains[source]
        for index, lines in enumerate(residues.values(), start=1):
            for line in lines:
                out.append(
                    f"{line[:6]}{serial:5d}{line[11:21]}{destination}"
                    f"{index:4d}{line[26]}{line[27:].rstrip()}\n"
                )
                serial += 1
        out.append(f"TER   {serial:5d}\n")
        serial += 1
        print(f"  chain {source} -> {destination}: {len(residues)} residues, "
              f"{sequence(residues)[:12]}{'...' if len(residues) > 12 else ''}")
    out.append("END\n")

    # The whole contract of this script: no atom moved. Compare the coordinate
    # columns as multisets -- order changes, values must not.
    src_atoms = [ln for residues in chains.values()
                 for lines in residues.values() for ln in lines]
    written = [ln for ln in out if ln.startswith("ATOM")]
    if len(written) != len(src_atoms):
        sys.exit(f"atom count changed: {len(src_atoms)} -> {len(written)}")
    if coord_set(written) != coord_set(src_atoms):
        sys.exit("COORDINATES CHANGED -- refusing to write. This script must "
                 "only relabel chains.")
    print(f"  {len(written)} atoms, coordinates identical")

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(out))
    print(f"wrote {dst}  (A=binder B=MHC C=peptide)")


if __name__ == "__main__":
    main()
