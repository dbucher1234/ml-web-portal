#!/usr/bin/env python
"""
Flask web server for MDCK Papp prediction using Open-source pipeline

Endpoints:
  GET  /            → renders index.html form
  POST /predict    → JSON results for SMILES or uploaded SDF
  POST /calculate_image → returns molecule image as base64
"""
import os
import base64
from io import BytesIO
import joblib
import pandas as pd
from flask import Flask, request, render_template, jsonify
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from utils_sasa import psa_from_mol

# Compute descriptors (same as train_model)
def compute_descriptors(mol: Chem.Mol) -> dict:
    mol = Chem.AddHs(mol)
    return {
        'MolWt': Chem.Descriptors.MolWt(mol),
        'LogP': Chem.Descriptors.MolLogP(mol),
        'TPSA': Chem.Descriptors.TPSA(mol),
        'NumHDonors': Chem.Descriptors.NumHDonors(mol),
        'NumHAcceptors': Chem.Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Chem.Descriptors.NumRotatableBonds(mol),
        'HeavyAtomCount': Chem.Descriptors.HeavyAtomCount(mol),
        'AromaticRingCount': Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
        'FormalCharge': Chem.rdmolops.GetFormalCharge(mol)
    }

# Load trained model pipeline once
global_pipeline = joblib.load(os.path.join('outputs', 'mdck_model.pkl'))
app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    # Accept SMILES or file upload
    smiles = request.form.get('smiles', '')
    sdf_file = request.files.get('sdf')
    try:
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            # generate 3D conformer
            params = AllChem.ETKDGv3()
            AllChem.EmbedMolecule(mol, params)
        elif sdf_file:
            suppl = Chem.SDMolSupplier(sdf_file.stream, removeHs=False)
            mol = suppl[0]
        else:
            return jsonify({'error': 'No input provided'}), 400

        # compute best conformer id (assumes single conf)
        conf_id = 0
        # compute 3D PSA
        psa3d = psa_from_mol(mol, conf_id)
        # compute descriptors
        desc = compute_descriptors(mol)
        # assemble feature vector
        data = {'3D_PSA': psa3d}
        data.update(desc)
        df = pd.DataFrame([data])
        # predict
        pred = global_pipeline.predict(df.values)[0]
        result = {
            '3D_PSA': round(psa3d, 2),
            **{k: round(v, 2) for k, v in desc.items()},
            'MDCK_logPapp': round(pred, 2)
        }
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/calculate_image', methods=['POST'])
def calculate_image():
    smiles = request.form.get('smiles', '')
    if not smiles:
        return jsonify({'error': 'No SMILES provided'}), 400
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({'error': 'Invalid SMILES'}), 400
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = BytesIO()
    img.save(buf, format='PNG')
    img_b64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    return jsonify({'mol_image': f"data:image/png;base64,{img_b64}"})

if __name__ == '__main__':
    app.run(debug=True)

