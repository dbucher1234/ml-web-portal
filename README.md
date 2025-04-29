# 🧪 3D PSA → Permeability → Web Portal

A three-step tutorial that starts with a **physics-based 3D polar surface area (PSA)** calculation, layers on a **machine-learning permeability model**, and ends with a **chemist-friendly Flask web app**.

---

## Why 3D PSA?

3D PSA measures the polar part of a molecule’s **solvent-accessible surface**.  
Unlike 2D TPSA, it captures shielding, intramolecular H-bonds and folding—key drivers of permeability in **beyond-Lipinski** space (macrocycles, PROTACs, peptides). A fast 3D estimate is therefore an excellent early PK filter.

---

## 🗺 Roadmap

| Step | Script / Folder | What it teaches | Speed |
|------|-----------------|-----------------|-------|
| **1. Compute 3D PSA** | `1_compute_psa.py` | Two routes:<br>• **QikProp** (Schrödinger) ⇒ reference 3D PSA<br>• **Open-source** (RDKit ETKDG + UFF → SASA triangulation) ⇒ ~10 Å² RMS vs QikProp | QP: 2–5 s / mol<br>OS: ≤0.5 s / mol |
| **2. Train ML model** | `2_train_mdck_model.py` | Gradient-Boost regressor that predicts **MDCK permeability** from 3D PSA (+ cLogP & MW) using a set of 328 cyclic peptides. | ~10 s total |
| **3. Build web portal** | `app.py`, `templates/`, `static/` | Flask app: paste a SMILES → server returns predicted 3D PSA and MDCK Papp in milliseconds. | ~50 ms / mol |

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

- Lawrenz M., Svensson M., Kato M. et al.  
  “A Computational Physics-based Approach to Predict Unbound Brain-to-Plasma Partition Coefficient, Kp,uu.”  
  *J. Chem. Inf. Model.* **2023**, 63, 12, 3786–3798.

- Möbitz H.  
  “Design Principles for Balancing Lipophilicity and Permeability in beyond Rule-of-5 Space.”  
  *ChemMedChem* **2023**, 18, e202300395.


