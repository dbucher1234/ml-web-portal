# 🧪 3D PSA → Permeability → Web Portal

Permeability is critical for large or beyond‐Lipinski compounds (macrocycles, peptides, PROTACs).  Their high polarity and size often hinder passive cellular uptake.

This project demonstrates a **fast 3D PSA calculation** + **machine‐learning MDCK permeability model**, culminating in a **chemist‐friendly prediction script** or Flask web portal.

---

## 📂 File Structure & Folders

```
ml-web-portal/
├── data/                   # Example input files
│   ├── penta_ala.smi       # SMILES for cyclic penta‑alanine
│   └── cycpep_training.csv # 40‑compound training set (SMILES, MDCK, 3D PSA)
│
├── outputs/                # Generated results and models
│   ├── penta_ala.sdf       # Lowest‑energy conformer 3D SDF (gen_conf.py)
│   ├── mdck_model.pkl      # Trained Random Forest model (train_model.py)
│   └── penta_ala_pred.csv  # Predicted MDCK values (predict_mdck.py)
│
├── src/                    # Python scripts for the workflow
│   ├── gen_conf.py         # Generate 3D conformers from SMILES
│   ├── utils_sasa.py       # Compute 3D PSA via FreeSASA on 3D SDF
│   ├── train_model.py      # Train RF regression on MDCK Papp
│   └── predict_mdck.py     # Predict MDCK Papp on new 3D SDF
│
└── environment.yml         # Conda environment with RDKit, FreeSASA, sklearn, etc.
```

---

## 🚀 Quick Tutorial: SMILES or SDF → MDCK Prediction

Below is a step-by-step guide to run the full pipeline on **any ligand** (SMILES or SDF).

1. **Create a Conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate ml_web
   ```

2. **Generate or provide a 3D conformer**

   If you have a SMILES:

   ```bash
   python src/gen_conf.py path/to/ligand.smi --output outputs/ligand.sdf
   ```

   If you already have a 3D SDF:

   ```bash
   cp path/to/ligand.sdf outputs/ligand.sdf
   ```

3. **Compute 3D PSA**

   ```bash
   python src/utils_sasa.py outputs/ligand.sdf
   ```

4. **Train the MDCK permeability model**

   ```bash
   python src/train_model.py \
       --input data/cycpep_training.csv \
       --output outputs/mdck_model.pkl
   ```

5. **Predict MDCK for your ligand**

   ```bash
   python src/predict_mdck.py \
       --model outputs/mdck_model.pkl \
       --input-sdf outputs/ligand.sdf \
       --output-csv outputs/ligand_pred.csv
   ```

---

## 🧪 Toy System Example: cyclic penta‑alanine

<p align="center">
  <img src="images/penta_ala.png" alt="Cyclic penta-alanine structure" />
</p>

SMILES (head‐to‐tail cyclized Ala₅):

```
O=C1[C@H](NC([C@H](NC([C@H](NC([C@H](NC([C@H](N1)C)=O)C)=O)C)=O)C)=O)C
```

**Example run outputs:**

1. **Generate 3D conformer**

```bash
python src/gen_conf.py data/penta_ala.smi --output outputs/penta_ala.sdf
# penta_ala: best conformer = 15, energy = 25.50 kcal/mol
```

2. **Compute 3D PSA**

```bash
python src/utils_sasa.py outputs/penta_ala.sdf
# penta_ala    PSA = 159.9 Å²
```

3. **Predict MDCK**

```bash
python src/predict_mdck.py \
    --model outputs/mdck_model.pkl \
    --input-sdf outputs/penta_ala.sdf
# penta_ala    Predicted_MDCK = -5.250217
```

> **Note:**\*\* MDCK Papp is reported in **log scale**. A log Papp (A→B) of –5 corresponds to 10⁻⁵ cm/s, indicating moderate cell permeability.

## 📈 Model Details & Performance

* **Model:** Random Forest regressor with hyperparameter grid over `n_estimators=[50,100]`, `max_depth=[5,10]` via 5‑fold CV.
* **Features:** 3D PSA (FreeSASA), MolWt, LogP, TPSA, H‑bond donors/acceptors, rotatable bonds, heavy atom count, aromatic ring count, formal charge.
* **Performance:** Achieves CV R² ≈ 0.41 and MAE ≈ 0.44 log units on the 40‑compound training set.

> Note: With only 40 samples, the model is proof‑of‑concept. Expanding the dataset can further improve accuracy.

---

## 🖥️ Web Portal Setup

A Flask-based web interface allows chemists to interact with the pipeline without writing code:

1. **Ensure the Flask app and templates** are in `src/` alongside `app.py`, along with `templates/` and `static/` directories.  The key file is `src/app.py`.
2. **Install dependencies** (RDKit, FreeSASA, Flask, etc.) via:

   ```bash
   conda activate ml_web
   pip install flask joblib pandas
   # or ensure environment.yml has flask and run `conda env update -f environment.yml`
   ```
3. **Run the server**:

   ```bash
   # from project root
   export FLASK_APP=src/app.py
   flask run
   # or simply
   python src/app.py
   ```
4. **Open your browser** at `http://localhost:5000`.
5. **Use the web form** to enter a SMILES string or upload a 3D SDF:

   * The app will generate (or accept) a 3D conformer
   * Compute 3D PSA and other descriptors
   * Return the **MDCK Papp** prediction in JSON or render on the page

This web portal leverages the same underlying scripts (`gen_conf.py`, `utils_sasa.py`, `predict_mdck.py`) for a seamless user experience.

## 📚 References

* Möbitz H. “Design Principles for Balancing Lipophilicity and Permeability in beyond Rule‑of‑5 Space.” *ChemMedChem* **2023**, 18, e202300395.
* Lawrenz M. et al. “A Computational Physics‑based Approach to Predict Unbound Brain‑to‑Plasma Partition Coefficient, Kp,uu.” *J. Chem. Inf. Model.* **2023**, 63(12), 3786–3798.

