# ML model → Web portal

Paste or draw any molecule (ChemDraw → SMILES)  
→ server builds a 3-D conformer  
→ ML model predicts its 3-D polar surface area (PSA)  
→ instant result + 2-D depiction in the browser.

[![Open in your browser](static/screenshot.png)](static/screenshot.png)

---

## 🚀 Quick start (local)

```bash
# 1. clone the repo
git clone https://github.com/dbucher1234/openfe-fep-aws.git
cd openfe-fep-aws/ml-web-portal

# 2. create the conda env
conda env create -f environment.yml
conda activate ml_web

# 3. launch
python app.py
# → http://127.0.0.1:5000

