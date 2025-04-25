import logging
from flask import Flask, request, render_template, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os
import subprocess
import pandas as pd
import joblib
import time
from io import BytesIO
import base64

app = Flask(__name__)

# Configure logging
logging.basicConfig(level=logging.DEBUG)

def generate_conformations(smiles, num_conformations=10):
    app.logger.debug(f"Generating conformations for SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)  # Add hydrogens

    params = AllChem.ETKDGv3()
    params.numThreads = 0  # Use all available threads
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations, params=params)

    return mol, ids

def optimize_conformations(mol, ids):
    app.logger.debug("Optimizing conformations")
    energies = []
    for conf_id in ids:
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
        energy = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id).CalcEnergy()
        energies.append((conf_id, energy))

    energies.sort(key=lambda x: x[1])
    lowest_energy_id = energies[0][0]

    return mol, lowest_energy_id

@app.route('/calculate_image', methods=['POST'])
def calculate_image():
    smiles = request.form['smiles']
    app.logger.debug(f"Received SMILES for image: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    mol_image = None
    if mol:
        mol_image = Draw.MolToImage(mol, size=(800, 800))  # Increased the image size two-fold
        buffered = BytesIO()
        mol_image.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        img_str = f"data:image/png;base64,{img_str}"
    return jsonify({'mol_image': img_str})

def save_conformation_to_sdf(mol, conf_id, output_file):
    app.logger.debug(f"Saving conformation to SDF: {output_file}")
    writer = Chem.SDWriter(output_file)
    writer.write(mol, confId=conf_id)
    writer.close()

def create_ligprep_inp():
    app.logger.debug("Creating ligprep.inp file")
    with open('ligprep.inp', 'w') as f:
        f.write("""
INPUT_FILE_NAME   output.sdf
MAX_ATOMS   500
FORCE_FIELD   16
USE_DESALTER   no
GENERATE_TAUTOMERS   no
DETERMINE_CHIRALITIES   no
IGNORE_CHIRALITIES   no
NUM_STEREOISOMERS   1
OUT_MAE   output.maegz
""")

def run_ligprep():
    app.logger.debug("Running LigPrep")
    schrodinger_path = '/opt/schrodinger2024-2'  # Ensure this path is correct
    
    # Log current directory contents
    app.logger.debug(f"Current directory: {os.getcwd()}")
    app.logger.debug(f"Files in current directory: {os.listdir('.')}")
    
    create_ligprep_inp()
    
    if not os.path.exists("output.sdf"):
        raise FileNotFoundError("The file output.sdf was not found before running LigPrep.")
    
    # Log environment variables
    app.logger.debug(f"Environment variables: {os.environ}")
    app.logger.debug(f"PATH: {os.environ['PATH']}")
    
    # Run LigPrep with Epik and capture output
    with open('ligprep_stdout.log', 'w') as stdout_log, open('ligprep_stderr.log', 'w') as stderr_log:
        try:
            result = subprocess.run([f"{schrodinger_path}/ligprep", "-inp", "ligprep.inp", "-epik"], check=True, stdout=stdout_log, stderr=stderr_log, text=True)
            app.logger.debug(f"LigPrep completed successfully.")
        except subprocess.CalledProcessError as e:
            app.logger.error(f"LigPrep failed with error: {e.stderr}")
            raise

def run_macrocycle_sampling(output_file):
    app.logger.debug("Running Macrocycle Sampling")
    schrodinger_path = '/opt/schrodinger2024-2'  # Ensure this path is correct
    macrocycle_command = (
        f"{schrodinger_path}/run macrocycle_conformational_sampling.py {output_file} -ffld S-OPLS -s WATER "
        f"-energy_window 10.0 -rmsd_cutoff 0.75 -sim_iterations 50 -iterations 50 -eigen_recalc global_min "
        f"-planar_torsion_sampling intermediate -j macrocycle_sampling -HOST localhost:24"
    )
    subprocess.run(macrocycle_command, shell=True, check=True)

def run_qikprop():
    app.logger.debug("Running QikProp")
    schrodinger_path = '/opt/schrodinger2024-2'  # Ensure this path is correct
    subprocess.run([f"{schrodinger_path}/qikprop", "output.maegz"], check=True)

def wait_for_file(filename, timeout=300):
    app.logger.debug(f"Waiting for file: {filename}")
    start_time = time.time()
    while not os.path.isfile(filename):
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Timeout waiting for file {filename}")
        time.sleep(5)

def extract_properties_from_csv(csv_file):
    app.logger.debug(f"Extracting properties from CSV: {csv_file}")
    df = pd.read_csv(csv_file)
    if df.empty:
        raise ValueError("The CSV file is empty.")
    
    properties = df[['mol_MW', 'QPlogPo/w', 'QPlogS', 'QPPCaco', 'QPPMDCK', 'PercentHumanOralAbsorption', 'PSA']].iloc[0].to_dict()
    properties['QPPMDCK'] /= 10
    properties['QPPCaco'] /= 10

    # Round properties to 2 significant digits
    properties = {k: round(v, 2) for k, v in properties.items()}

    return properties

def prepare_data_for_model(csv_file):
    app.logger.debug(f"Preparing data for model from CSV: {csv_file}")
    df = pd.read_csv(csv_file)
    if df.empty:
        raise ValueError("The CSV file is empty.")
    
    data_prepared = df[['QPlogPo/w', 'PSA', 'mol_MW']]
    data_prepared.columns = ['qp', 'psa', 'MW']
    return data_prepared

def calculate_properties(smiles):
    app.logger.debug("Starting calculation of properties")
    mol, ids = generate_conformations(smiles)
    mol, lowest_energy_id = optimize_conformations(mol, ids)
    save_conformation_to_sdf(mol, lowest_energy_id, "output.sdf")

    run_ligprep()
    wait_for_file("output.maegz")
    run_macrocycle_sampling("output.sdf")
    run_qikprop()
    
    wait_for_file("output.CSV")
    properties = extract_properties_from_csv("output.CSV")
    
    return properties

def predict_mdck():
    app.logger.debug("Predicting MDCK")
    csv_file = 'output.CSV'
    data_prepared = prepare_data_for_model(csv_file)
    model = joblib.load('trained_model.pkl')
    prediction = model.predict(data_prepared)[0]
    # Rescale predictions
    rescaled_prediction = (prediction + 5.261) / 0.1479
    return 10 ** rescaled_prediction * 1000000

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    smiles = request.form['smiles']
    app.logger.debug(f"Received SMILES for calculation: {smiles}")
    
    try:
        # Calculate properties
        properties = calculate_properties(smiles)
        
        # Predict using the trained model
        model_prediction = predict_mdck()
        model_prediction_rounded = round(model_prediction, 2)
        
        # Combine results
        results = {
            **properties,
            'MDCK_Prediction': model_prediction_rounded
        }
        
        return jsonify(results)
    except Exception as e:
        app.logger.error(f"Error during calculation: {e}")
        return jsonify({'error': str(e)})

if __name__ == '__main__':
    app.run(debug=True)

