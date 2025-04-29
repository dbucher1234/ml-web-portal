#!/usr/bin/env python
"""
gen_conf.py
===========

SMILES → 10 low-energy conformers using RDKit ETKDG + UFF.
Writes a multi-conformer SDF for each ligand.

Usage
-----
python gen_conf.py ligands.smi  --out conformers/
"""

import argparse, os
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

def build_conformers(smiles: str, n: int = 10):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, numConfs=n, params=AllChem.ETKDGv3())
    for cid in range(mol.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(mol, confId=cid)
    return mol

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("smifile", help="SMILES file (one per line)")
    parser.add_argument("--out", default="conformers",
                        help="Output folder for .sdf files")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    with open(args.smifile) as fh:
        for line in fh:
            smi, *name = line.strip().split()
            name = name[0] if name else f"lig_{fh.tell()}"
            mol = build_conformers(smi)
            path = os.path.join(args.out, f"{name}.sdf")
            w = SDWriter(path)
            for cid in range(mol.GetNumConformers()):
                mol.SetProp("_Name", f"{name}_conf{cid}")
                w.write(mol, confId=cid)
            w.close()
            print(f"Wrote {path}")

