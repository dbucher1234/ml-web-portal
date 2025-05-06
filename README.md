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

## 🚀 Quick Tutorial: SMILES → MDCK Prediction

Below is a step‑by‑step guide to run the full pipeline on **penta‑ala**.

1. **Create a Conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate ml_web
   ```

2. **Generate a 3D conformer**

   ```bash
   python src/gen_conf.py data/penta_ala.smi --output outputs/penta_ala.sdf
   ```

   * Uses RDKit ETKDG + UFF to sample 20 conformers
   * Selects the lowest‐energy structure

3. **Compute 3D PSA (polar surface area)**

   ```bash
   python src/utils_sasa.py outputs/penta_ala.sdf
   ```

   * Calls FreeSASA to sum solvent‐accessible surface on N/O atoms

4. **Train the MDCK permeability model**

   ```bash
   python src/train_model.py \
       --input data/cycpep_training.csv \
       --output outputs/mdck_model.pkl
   ```

   * Loads training CSV (SMILES, MDCK, 3D PSA)
   * Computes eight descriptors (MolWt, LogP, TPSA, H‑bond donors/acceptors, rotatable bonds, heavy atom count, aromatic rings, formal charge)
   * Trains a Random Forest (5‑fold CV) and reports CV R² & MAE
   * Saves the trained pipeline to `outputs/mdck_model.pkl`

5. **Predict MDCK on a new peptide**

   ```bash
   python src/predict_mdck.py \
       --model outputs/mdck_model.pkl \
       --input-sdf outputs/penta_ala.sdf \
       --output-csv outputs/penta_ala_pred.csv
   ```

   * Loads your 3D SDF, recomputes PSA + descriptors
   * Outputs predicted MDCK Papp values to CSV

---

## 🧪 Toy System Example: cyclic penta‑alanine

SMILES (head‐to‐tail cyclized Ala₅):

```
O=C1[C@H](NC([C@H](NC([C@H](NC([C@H](NC([C@H](N1)C)=O)C)=O)C)=O)C)=O)C
```

This drives the steps above and yields a predicted MDCK log Papp.

---

## 📈 Model Details & Performance

* **Model:** Random Forest regressor with hyperparameter grid over `n_estimators=[50,100]`, `max_depth=[5,10]` via 5‑fold CV.
* **Features:** 3D PSA (FreeSASA), MolWt, LogP, TPSA, H‑bond donors/acceptors, rotatable bonds, heavy atom count, aromatic ring count, formal charge.
* **Performance:** Achieves CV R² ≈ 0.41 and MAE ≈ 0.44 log units on the 40‑compound training set.

> Note: With only 40 samples, the model is proof‑of‑concept. Expanding the dataset or exploring gradient boosting and kernel methods can further improve accuracy.

---

## 📚 References

* Möbitz H. “Design Principles for Balancing Lipophilicity and Permeability in beyond Rule‑of‑5 Space.” *ChemMedChem* **2023**, 18, e202300395.
* Lawrenz M. et al. “A Computational Physics‑based Approach to Predict Unbound Brain‑to‑Plasma Partition Coefficient, Kp,uu.” *J. Chem. Inf. Model.* **2023**, 63(12), 3786–3798.

