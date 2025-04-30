python - <<'PY'
from rdkit import Chem
import sys, os
sdf = sys.argv[1] if len(sys.argv) > 1 else "test_conf/acetone.sdf"
mol = Chem.SDMolSupplier(sdf, removeHs=False)[0]
print("Molecule read:", bool(mol))
print("Conformers:", mol.GetNumConformers() if mol else 0)
if mol and mol.GetNumConformers():
    conf = mol.GetConformer()
    print("First atom coords:", conf.GetAtomPosition(0))
PY  test_conf/acetone.sdf

