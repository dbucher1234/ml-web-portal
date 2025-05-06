#!/usr/bin/env python
"""
predict_mdck.py

Use a trained Random Forest model to predict MDCK Papp for molecules with 3D coordinates.

Usage:
    python src/predict_mdck.py \
        --model ../outputs/mdck_model.pkl \
        --input-sdf ../outputs/penta_ala.sdf \
        [--output-csv predictions.csv]

Outputs predictions to console or CSV.
"""
import argparse
import joblib
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import numpy as np
from utils_sasa import psa_from_mol


def compute_descriptors(mol: Chem.Mol) -> dict:
    """
    Compute property-based descriptors for an RDKit Mol with explicit Hs.
    """
    return {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'HeavyAtomCount': Descriptors.HeavyAtomCount(mol),
        'AromaticRingCount': rdMolDescriptors.CalcNumAromaticRings(mol),
        'FormalCharge': Chem.rdmolops.GetFormalCharge(mol)
    }


def load_molecules(sdf_path: Path) -> list[Chem.Mol]:
    """Load molecules (with explicit Hs and 3D coords) from an SDF file."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    return mols


def main(model_path: Path, input_sdf: Path, output_csv: Path = None) -> None:
    # Load trained pipeline
    pipeline = joblib.load(model_path)

    # Read molecules
    mols = load_molecules(input_sdf)
    if not mols:
        print(f"No valid molecules found in {input_sdf}")
        return

    # Prepare data
    records = []
    for i, mol in enumerate(mols):
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol_{i}'
        # 3D PSA
        psa3d = psa_from_mol(mol, conf_id=0)
        # property descriptors
        props = compute_descriptors(mol)
        # assemble feature vector in same order as training
        data = {'3D_PSA': psa3d}
        data.update(props)
        records.append((name, data))

    # Create DataFrame
    names = [r[0] for r in records]
    feat_df = pd.DataFrame([r[1] for r in records])

    # Predict
    preds = pipeline.predict(feat_df.values)

    # Output
    out_df = pd.DataFrame({'Name': names, 'Predicted_MDCK': preds})
    if output_csv:
        out_df.to_csv(output_csv, index=False)
        print(f"Predictions written to {output_csv}")
    else:
        print(out_df.to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Predict MDCK Papp using trained model and SDF")
    parser.add_argument('--model', type=Path,
                        default=Path(__file__).parents[1] / 'outputs' / 'mdck_model.pkl',
                        help='Path to trained model pickle')
    parser.add_argument('--input-sdf', type=Path,
                        default=Path(__file__).parents[1] / 'outputs' / 'penta_ala.sdf',
                        help='Input SDF file with 3D coordinates')
    parser.add_argument('--output-csv', type=Path,
                        help='Optional path to write CSV of predictions')
    args = parser.parse_args()
    main(args.model, args.input_sdf, args.output_csv)

