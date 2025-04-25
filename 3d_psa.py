import os
import subprocess
import numpy as np
import time
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

def generate_multiple_conformations(smiles, num_conformations=10):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations, params=params)
    for conf_id in conf_ids:
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
    return mol, conf_ids

def export_conformations_to_sdf(mol, conf_ids, filename):
    writer = SDWriter(filename)
    for conf_id in conf_ids:
        writer.write(mol, confId=conf_id)
    writer.close()

def run_qikprop(input_file):
    schrodinger_path = "/opt/schrodinger2024-2"
    qikprop_command = f"{schrodinger_path}/qikprop {input_file}"
    
    print(f"Running QikProp with command: {qikprop_command}")
    result = subprocess.run(qikprop_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    print("QikProp stdout:")
    print(result.stdout.decode())
    print("QikProp stderr:")
    print(result.stderr.decode())
    
    # List directory contents to check generated files
    print("Directory contents after QikProp run:")
    for f in os.listdir('.'):
        print(f)
    
    # Wait for the output file to be created
    output_file_base = os.path.splitext(input_file)[0]
    output_file = f"{output_file_base}.qpsa"
    
    wait_time = 0
    max_wait_time = 60  # Maximum wait time in seconds
    while not os.path.exists(output_file) and wait_time < max_wait_time:
        time.sleep(5)
        wait_time += 5
        print(f"Waiting for the output file: {output_file}. Elapsed time: {wait_time} seconds")
    
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"QikProp did not create the expected output file: {output_file}")
    
    return output_file

def parse_qpsa_output(output_file):
    psa_values = []
    with open(output_file, 'r') as file:
        for line in file:
            if "PSA" in line:
                psa_value = float(line.split()[-1])
                psa_values.append(psa_value)
    return psa_values

def calculate_statistics(psa_values):
    average_psa = np.mean(psa_values)
    std_dev_psa = np.std(psa_values)
    return average_psa, std_dev_psa

def main():
    # Ask for SMILES input
    smiles = input("Enter the SMILES string: ")
    sdf_output = "molecule_conformations.sdf"
    num_conformations = 10

    # Generate multiple 3D conformations
    mol, conf_ids = generate_multiple_conformations(smiles, num_conformations)
    
    # Export conformations to SDF
    export_conformations_to_sdf(mol, conf_ids, sdf_output)
    
    # Run QikProp
    qikprop_output = run_qikprop(sdf_output)
    
    # Parse and print the PSA results
    psa_values = parse_qpsa_output(qikprop_output)
    average_psa, std_dev_psa = calculate_statistics(psa_values)
    print(f"3D PSA Average: {average_psa:.2f}")
    print(f"3D PSA Standard Deviation: {std_dev_psa:.2f}")

if __name__ == "__main__":
    main()

