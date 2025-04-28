# ML model → Web portal

Why 3-D PSA matters?

Polar surface area (PSA) is the part of a molecule that can hydrogen-bond with water. A 3-D PSA uses the molecule’s real 3-D shape, so it captures shielding, intramolecular H-bonds, and folding—factors that strongly influence permeability and oral exposure, especially for non-Lipinski “beyond-Rule-of-5” compounds. Classic 2-D TPSA is quick but ignores conformation, so it often under- or over-estimates permeability for flexible or macrocyclic molecules. A fast, accurate 3-D PSA prediction is therefore a handy early PK filter.

What does this demo do?
Give it a SMILES and you get a 3-D polar surface-area (PSA) estimate two ways:

1) Full 3-D route (optional) – builds real 3-D conformers and, if a Schrödinger install is detected, runs QikProp for the “true” PSA. Accurate but a few seconds per molecule.
  
2) Fast ML fallback – computes cheap 2-D TPSA, then adds a machine-learned Δ-correction that was trained on thousands of QikProp results. No 3-D step, returns a PSA within ≈5 Å² of QikProp in milliseconds.
Use the first when you need highest fidelity; use the second for instant, chemist-friendly screening.

# Outline of the Tool

Paste or draw any molecule (ChemDraw → SMILES)  
→ server builds a 3-D conformer  
→ ML model predicts its 3-D polar surface area (PSA)  
→ instant result + 2-D depiction in the browser.

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

