#!/usr/bin/env python
"""
3d_psa.py
=========

Educational utility functions for the **ML-powered Web Portal** demo.

Goal
----
1. Take an input SMILES (from ChemDraw or text box).
2. Generate multiple 3-D conformations with RDKit (ETKDG + UFF).
3. Optionally run QikProp (if SCHRODINGER env var is set) to compute
   “true” 3-D polar surface area (PSA).
4. Load a pretrained machine-learning model (joblib) that predicts
   3-D PSA directly from molecular descriptors (e.g., 2-D TPSA + Δ-learned
   correction).
5. Return both *predicted* and *reference* PSA values for the web app.

Notes
-----
* The script is **stand-alone**; import its functions in `app.py`.
* Works fine with open-source RDKit only; QikProp call is skipped if
  Schrödinger tools are not available.
"""

import os
import subprocess
import numpy as np
from typing import Tuple, List

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

# ---------------------------------------------------------------------
# 1 ▸ 3-D conformer generation
# ---------------------------------------------------------------------
def generate_multiple_conformations(
    smiles: str,
    num_conformations: int = 10
) -> Tuple[Chem.Mol, List[int]]:
    """
    Embed `num_conformations` 3-D structures for a molecule.

    Parameters
    ----------
    smiles : str
        Input SMILES string.
    num_conformations : int, optional
        How many ETKDG conformers to generate (default = 10).

    Returns
    -------
    mol : rdkit.Chem.Mol
        RDKit molecule with hydrogens and 3-D coordinates.
    conf_ids : list[int]
        RDKit conformer IDs that were successfully embedded.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations,
                                          params=params)
    # Local geometry optimisation (UFF) for each conformer
    for cid in conf_ids:
        AllChem.UFFOptimizeMolecule(mol, confId=cid)
    return mol, conf_ids


# ---------------------------------------------------------------------
# 2 ▸ Write conformers to SDF (helper for QikProp or visualisation)
# ---------------------------------------------------------------------
def write_conformers_to_sdf(mol: Chem.Mol, path: str) -> None:
    """Save all conformers in a single SDF file."""
    writer = SDWriter(path)
    for cid in range(mol.GetNumConformers()):
        mol.SetProp("_Name", f"conf_{cid}")
        writer.write(mol, confId=cid)
    writer.close()


# ---------------------------------------------------------------------
# 3 ▸ Optional: call QikProp for reference 3-D PSA
# ---------------------------------------------------------------------
def run_qikprop(input_sdf: str) -> str:
    """
    Run Schrödinger QikProp on the provided SDF (if available).

    Returns
    -------
    str : path to the generated `.qpsa` output file.

    Raises
    ------
    RuntimeError if QikProp is not found or fails to execute.
    """
    schrodinger_path = os.environ.get("SCHRODINGER")
    if not schrodinger_path:
        raise RuntimeError("SCHRODINGER env var not set; skipping QikProp")

    cmd = [os.path.join(schrodinger_path, "qikprop"), input_sdf]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"QikProp failed:\n{result.stderr.decode()}")

    base = os.path.splitext(input_sdf)[0]
    return f"{base}.qpsa"


# ---------------------------------------------------------------------
# 4 ▸ Convenience wrapper: SMILES → conformers → (optional) QikProp
# ---------------------------------------------------------------------
def prepare_ligand(smiles: str, use_qikprop: bool = False) -> Tuple[str, str]:
    """
    Full pipeline for a single ligand.

    Returns
    -------
    Tuple[input_sdf, qpsa_path_or_none]
    """
    mol, _ = generate_multiple_conformations(smiles)
    sdf_path = "current_ligand.sdf"
    write_conformers_to_sdf(mol, sdf_path)

    if use_qikprop:
        try:
            qpsa_file = run_qikprop(sdf_path)
        except RuntimeError as err:
            print(f"[WARN] {err}")
            qpsa_file = ""
    else:
        qpsa_file = ""

    return sdf_path, qpsa_file


# ---------------------------------------------------------------------
# 5 ▸ Entry-point for command-line testing
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate conformers and (optionally) run QikProp."
    )
    parser.add_argument("smiles", help="SMILES string for the ligand")
    parser.add_argument("--qikprop", action="store_true",
                        help="Run QikProp if SCHRODINGER path is set")
    args = parser.parse_args()

    sdf, qpsa = prepare_ligand(args.smiles, use_qikprop=args.qikprop)
    print(f"SDF written to: {sdf}")
    if qpsa:
        print(f"QikProp output: {qpsa}")

