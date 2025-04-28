# Why predict 3D PSA?

Polar surface area (PSA) is the part of a molecule that interacts with water. A 3D PSA uses the molecule’s real 3D shape, so it captures shielding, intramolecular H-bonds, and folding—factors that strongly influence permeability and oral exposure, especially for non-Lipinski “beyond-Rule-of-5” compounds. Classic 2D TPSA is quick but ignores conformation, so it often under- or over-estimates permeability for flexible or macrocyclic molecules. A fast, accurate 3D PSA prediction is therefore a handy early PK filter.

What does this demo do?
Give it a SMILES and you get a 3D polar surface-area (PSA) estimate two ways:

1) Full 3D route (optional) – builds real 3D conformers and, if a Schrödinger install is detected, runs QikProp for the “true” PSA. Accurate but a few seconds per molecule.
  
2) Fast ML fallback – computes cheap 2D TPSA, then adds a machine-learned Δ-correction that was trained on thousands of QikProp results. No 3D step, returns a PSA within ≈5 Å² of QikProp in milliseconds.
Use the first when you need highest fidelity; use the second for instant, chemist-friendly screening.

# ML model → Web portal

Paste or draw any molecule (ChemDraw → SMILES)  
→ server builds a 3D conformer  
→ ML model predicts its 3D polar surface area (PSA)  
→ instant result + 2D depiction in the browser.

[![Open in your browser](static/screenshot.png)](static/screenshot.png)

---

## 🚀 Quick start (local)

```bash
# 1. clone the repo
git clone https://github.com/dbucher1234/ml-web-portal.git
cd ml-web-portal

# 2. create the conda env
conda env create -f environment.yml
conda activate ml_web

# 3. launch
python app.py
# → http://127.0.0.1:5000

