# 🧪 3D PSA → Permeability → Web Portal

This tutorial starts with a **physics-based 3D polar surface area (PSA)** calculation, layers on a **machine-learning permeability model**, and ends with a **chemist-friendly Flask web app**.

---

## Why 3D PSA?

Permeability becomes critical for compounds above ~650 Da or outside Lipinski space because their large size and high polarity hinder passive diffusion across cell membranes, making efficient cellular uptake the primary barrier to efficacy.

3D PSA measures the polar part of a molecule’s solvent-accessible surface. Unlike 2D TPSA, it captures **shielding, intramolecular H-bonds and folding—key drivers of permeability in beyond-Lipinski space** (e.g. macrocycles, PROTACs, peptides).

To address these permeability challenges, several groups have turned to physics-based 3D descriptors like 3D PSA and solvation energy (E-sol), each implementing distinct workflows to capture conformational and electronic effects. Möbitz et al. generated 3D PSA by first producing conformers with OpenEye Omega and clustering them with RDKit. Each representative conformer then underwent a single-point COSMO QM calculation, after which the polar surface area was obtained by summing the solvent-accessible surface where the atomic charge density exceeded |q| > 0.002 e Å⁻². These 3D PSA values were subsequently used to train a permeability model on 114 proprietary beyond-Rule-of-5 compounds.

Additionally, Lawrenz et al. employed a Schrödinger workflow, in which they built 3D structures with LigPrep and sampled conformations in MacroModel. For the lowest-energy conformers, they computed single-point QM calculations both in gas phase and with an implicit-solvent model. The difference between the two energies provided the 3D solvation energy, E-sol, which served as the key descriptor for a trained permeability model.

In this example, we propose a similar approach: 

---

## 🗺 Roadmap

| Step | Script / Folder | What it does | 
|------|-----------------|-----------------|
| **1. SMILES → conformers** | `gen_conf.py` | Use Open-source RDKit ETKDG to generate 10 low energy conformations per ligand.
| **2. Compute 3D PSA** | `compute_psa.py` | QikProp (Schrödinger) ⇒ 3D PSA<br> |
| **3. Train ML model** | `train_mdck_model.py` | Gradient-Boost regressor that predicts **MDCK permeability** from 3D PSA, in our example using a set of 328 cyclic peptides from http://cycpeptmpdb.com |
| **4. Build web portal** | `app.py`, `templates/`, `static/` | Flask app: paste a SMILES in server to return predicted 3D PSA and MDCK Papp in <1 second. |

---

## 🚀 Quick Start (local)

```bash
# clone
git clone https://github.com/dbucher1234/ml-web-portal.git
cd ml-web-portal

# create conda env
conda env create -f environment.yml
conda activate ml_web

# 1️⃣  compute PSA (open-source route)
python 1_compute_psa.py data/ligands.smi --method open

# 2️⃣  train permeability model
python 2_train_mdck_model.py

# 3️⃣  launch web portal
python app.py           # → http://127.0.0.1:5000

```
---

## 📚 References

- Möbitz H.  
  “Design Principles for Balancing Lipophilicity and Permeability in beyond Rule-of-5 Space.”  
  *ChemMedChem* **2023**, 18, e202300395.
  
- Lawrenz M., Svensson M., Kato M. et al.  
  “A Computational Physics-based Approach to Predict Unbound Brain-to-Plasma Partition Coefficient, Kp,uu.”  
  *J. Chem. Inf. Model.* **2023**, 63, 12, 3786–3798.




