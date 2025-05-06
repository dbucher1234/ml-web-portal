#!/usr/bin/env python
"""
utils_sasa.py – quiet 3-D polar surface area (Å²) calculator
-----------------------------------------------------------
• Input  : SDF or PDB with 3-D coordinates and explicit hydrogens
• Output : one line per conformer with title and PSA value
• Polar atoms counted: N and O (edit POLAR_ELEMS to include others)
"""

from __future__ import annotations
import argparse, os, pathlib, sys, tempfile
from rdkit import Chem
import freesasa
import contextlib

# Configuration
PROBE_RADIUS = 1.4            # Å (match QikProp default)
POLAR_ELEMS  = {"N", "O"}   # add "S", "Cl", etc. if needed

@contextlib.contextmanager
def silence_stderr():
    """Suppress C-level stderr for FreeSASA calls."""
    fd = os.open(os.devnull, os.O_RDWR)
    old_stderr = os.dup(2)
    try:
        os.dup2(fd, 2)
        yield
    finally:
        os.dup2(old_stderr, 2)
        os.close(fd)
        os.close(old_stderr)


def psa_from_mol(mol: Chem.Mol, conf_id: int = 0) -> float:
    """Compute 3-D polar surface area (Å²) for a given conformer."""
    pdb_block = Chem.MolToPDBBlock(mol, confId=conf_id)
    with tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=False) as tmp:
        tmp.write(pdb_block)
        tmp.flush()
        pdb_path = tmp.name

    try:
        with silence_stderr():
            structure = freesasa.Structure(
                pdb_path,
                options={"hydrogen": True, "hetatm": True}
            )
        with silence_stderr():
            params = freesasa.Parameters({"probe-radius": PROBE_RADIUS})
            result = freesasa.calc(structure, params)

        psa = 0.0
        for idx, atom in enumerate(mol.GetAtoms()):
            if atom.GetSymbol() in POLAR_ELEMS:
                psa += result.atomArea(idx)
        return psa

    finally:
        os.remove(pdb_path)


def load_mols(path: pathlib.Path) -> list[Chem.Mol]:
    ext = path.suffix.lower()
    if ext == ".sdf":
        return [
            m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m
        ]
    if ext == ".pdb":
        m = Chem.MolFromPDBFile(str(path), removeHs=False)
        return [m] if m else []
    raise ValueError(f"Unsupported file type: {ext}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quiet 3-D polar surface area (Å²) calculator"
    )
    parser.add_argument(
        "infile",
        help="Input .sdf or .pdb with 3-D coords and hydrogens"
    )
    args = parser.parse_args()

    path = pathlib.Path(args.infile).expanduser().resolve()
    mols = load_mols(path)
    if not mols:
        sys.exit(
            "ERROR: could not read molecule(s) or no conformer found."
        )

    for mol in mols:
        title = (
            mol.GetProp("_Name") if mol.HasProp("_Name")
            else path.stem
        )
        psa = psa_from_mol(mol, conf_id=0)
        print(f"{title}\tPSA = {psa:.1f} Å²")


if __name__ == "__main__":
    main()

