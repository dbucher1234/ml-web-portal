# 🧪 3D PSA → Permeability → Web Portal

Permeability becomes critical for compounds above ~650 Da or outside Lipinski space (e.g. macrocycles, PROTACs, peptides), because their large size and high polarity hinder passive diffusion across cell membranes, making efficient cellular uptake one of the primary barrier to efficacy.

This tutorial starts with a **physics-based 3D polar surface area (PSA)** calculation, layers on a **machine-learning permeability model**, and ends with a **chemist-friendly Flask web app**.

---

## Why 3D PSA?

3D PSA measures the polar part of a molecule’s solvent-accessible surface. Unlike 2D TPSA, it captures **shielding, intramolecular H-bonds and folding—key drivers of permeability in beyond-Lipinski space**.

To address these permeability challenges, several groups have turned to physics-based 3D descriptors like 3D PSA and solvation energy (E-sol), each implementing distinct workflows to capture conformational and electronic effects. Möbitz et al. generated 3D PSA by first producing conformers with OpenEye Omega and clustering them with RDKit. Each representative conformer then underwent a single-point COSMO QM calculation, after which the polar surface area was obtained by summing the solvent-accessible surface where the atomic charge density exceeded |q| > 0.002 e Å⁻². These 3D PSA values were subsequently used to train a permeability model on 114 proprietary beyond-Rule-of-5 compounds.

Additionally, Lawrenz et al. employed a Schrödinger workflow, in which they built 3D structures with LigPrep and sampled conformations in MacroModel. For the lowest-energy conformers, they computed single-point QM calculations both in gas phase and with an implicit-solvent model. The difference between the two energies provided the 3D solvation energy, E-sol, which served as the key descriptor for a trained permeability model.

Here, we propose a similar approach, optimized for speed: 

---

## 🗺 Roadmap

| Step | Script / Folder | What it does | 
|------|-----------------|-----------------|
| **1. SMILES → conformers**     | `gen_conf.py` | Open-source RDKit (ETKDG + UFF) is used to sample 20 conformers and select the lowest-energy one, providing a fast but replaceable method for generating reasonable 3D structures. 
| **2. Compute 3D PSA**          | `compute_psa.py` | QikProp (Schrödinger) ⇒ to compute 3D PSA<br> (RDKIT can be used as an Open-source alternative). |
| **3. Build ML model**          | `mdck_model.py` | In our example, an **MDCK permeability** model for cyclic peptides was trained with Schrodinger AutoQSAR, using 3D PSA and standard property-based descriptors. The training set was based on 328 cyclic peptides from http://cycpeptmpdb.com |
| **4. Build web portal**        | `app.py`, `templates/`, `static/` | Flask app: paste a SMILES in server to return predicted 3D PSA and MDCK Papp in <1 sec. |

---

# 🧪 Toy System: cyclic L-Ala₅

We will use this as a test-case. The SMILES for a cyclic penta-alanine (five alanines linked head-to-tail in a ring) is: <pre markdown="1"> ```O=C1[C@H](NC([C@H](NC([C@H](NC([C@H](NC([C@H](N1)C)=O)C)=O)C)=O)C)=O)C ``` </pre>

Alternatively, one can draw it in a sketcher like ChemDraw and right-click + copy as smile. 

<pre markdown="1"> ```python gen_conf.py Ala5.smi  ``` </pre>
penta_ala: best conformer = 13, energy = 25.52 kcal/mol

We will now calculate the 3D PSA for this conformer:
<pre markdown="1"> ```python 3d_psa.py ``` </pre>
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




