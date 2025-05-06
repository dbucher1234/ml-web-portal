#!/usr/bin/env python
"""
gen_conf.py  —  SMILES → lowest-energy conformer (SDF)

Generates 3D conformers using RDKit ETKDG, optimizes them with UFF,
and writes only the lowest-energy conformer to SDF.

Usage:
  python gen_conf.py ligands.smi --out output_folder

Dependencies:
  - RDKit
"""

import os
import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_lowest_energy_conformer(smiles: str, name: str, out_dir: Path, n_confs: int = 20):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    # Optimize and score all conformers
    energies = []
    for cid in conf_ids:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        energies.append((cid, energy))

    # Select lowest-energy conformer
    best_cid, best_energy = min(energies, key=lambda x: x[1])
    print(f"{name}: best conformer = {best_cid}, energy = {best_energy:.2f} kcal/mol")

    mol.SetProp("_Name", name)
    writer = Chem.SDWriter(str(out_dir / f"{name}.sdf"))
    writer.write(mol, confId=best_cid)
    writer.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smi", help="Input SMILES file (e.g. ligands.smi)")
    parser.add_argument("--out", default="conf_out", help="Output folder [default: conf_out]")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)

    with open(args.smi) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 1:
                smiles, name = parts[0], "mol"
            else:
                smiles, name = parts[0], parts[1]
            generate_lowest_energy_conformer(smiles, name, out_dir)


if __name__ == "__main__":
    main()

