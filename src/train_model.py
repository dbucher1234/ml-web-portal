#!/usr/bin/env python
"""
train_model.py

Train a Random Forest regression model for MDCK Papp using 3D PSA and extended property-based descriptors.
Usage:
    python src/train_model.py \
        --input ../data/cycpep_training.csv \
        --output ../outputs/mdck_model.pkl
"""
import argparse
import joblib
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict, GridSearchCV
from sklearn.metrics import r2_score, mean_absolute_error


def compute_descriptors(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
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


def build_pipeline() -> Pipeline:
    """Construct a pipeline with scaling and Random Forest with hyperparameter tuning."""
    steps = [('scale', StandardScaler())]
    rf = RandomForestRegressor(random_state=42)
    param_grid = {
        'n_estimators': [50, 100],
        'max_depth': [5, 10]
    }
    # Use 5-fold CV for grid search
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    grid = GridSearchCV(rf, param_grid, cv=kf, scoring='r2')
    steps.append(('model', grid))
    return Pipeline(steps)


def main(input_path: Path, output_path: Path) -> None:
    # Load data
    df = pd.read_csv(input_path)
    desc_list = [compute_descriptors(smi) for smi in df['SMILES']]
    desc_df = pd.DataFrame(desc_list)

    # Feature matrix: 3D PSA + descriptors
    X = pd.concat([df[['3D_PSA']].reset_index(drop=True), desc_df], axis=1).values
    y = df['MDCK'].values

    # Build and evaluate pipeline
    pipeline = build_pipeline()
    # 5-fold CV for performance estimation
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    y_pred = cross_val_predict(pipeline, X, y, cv=kf)
    r2 = r2_score(y, y_pred)
    mae = mean_absolute_error(y, y_pred)
    print(f"CV R² = {r2:.2f}")
    print(f"CV MAE = {mae:.2f} log units")

    # Fit on full data and inspect best params
    pipeline.fit(X, y)
    best = pipeline.named_steps['model']
    print(f"Best RF params: {best.best_params_}")

    # Save pipeline
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(pipeline, output_path)
    print(f"Model pipeline saved to {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Train MDCK Papp Random Forest model")
    parser.add_argument('--input', type=Path,
                        default=Path(__file__).parents[1] / 'data' / 'cycpep_training.csv',
                        help='Path to input CSV file')
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).parents[1] / 'outputs' / 'mdck_model.pkl',
                        help='Path to save trained pipeline')
    args = parser.parse_args()
    main(args.input, args.output)

