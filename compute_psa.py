#!/usr/bin/env python
"""
compute_psa.py
==============

Compute 3D polar surface area (PSA) for an SDF:

1. QikProp route  (if --method qikprop and $SCHRODINGER present)
2. Open-source route: load SDF, triangulate SASA with RDKit + FreeSASA
   and sum polar atoms' surface.

Usage
-----
python compute_psa.py conformers/lig_A.sdf      # open-source default
python compute_psa.py conformers/lig_A.sdf --method qikprop
"""

import os, subprocess, argparse
import freesasa, rdkit.Chem as Chem

def psa_opensource(sdf:str) -> float:
    mol = Chem.SDMolSupplier(sdf, removeHs=False)[0]
    struct = freesasa.Structure(sdf)
    result = freesasa.calc(struct)
    psa = 0.0
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() in (7, 8):                   # N or O
            psa += result.atomArea(i)
    return psa

def psa_qikprop(sdf:str) -> float:
    sch = os.environ.get("SCHRODINGER")
    if not sch: raise RuntimeError("SCHRODINGER env var not set")
    subprocess.run([os.path.join(sch,"qikprop"), sdf], check=True)
    qprop = os.path.splitext(sdf)[0] + ".qprop"
    with open(qprop) as fh:
        for line in fh:
            if line.startswith("  PSA"):
                return float(line.split()[2])
    raise ValueError("PSA not found in .qprop")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("sdf", help="Input SDF with conformers")
    p.add_argument("--method", choices=["open", "qikprop"],
                   default="open", help="PSA engine")
    args = p.parse_args()

    if args.method == "qikprop":
        psa = psa_qikprop(args.sdf)
    else:
        psa = psa_opensource(args.sdf)

    print(f"3D PSA = {psa:.1f} Å²  ({args.method})")

