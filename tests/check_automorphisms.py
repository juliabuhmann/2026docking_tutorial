"""Every ligand's symmetry matchings must be true graph automorphisms: include the identity, map
each atom to one of the same element, and map bonds onto bonds. Run after scripts/build_data.py."""
import json
import re
from pathlib import Path

from rdkit import Chem

DATA = Path(__file__).resolve().parents[1] / "docs" / "data"
failures = 0
for f in sorted(DATA.glob("prot_*.js")):
    d = json.loads(re.search(r"PW\.data\([^,]+,\s*[^,]+,\s*(.*)\);\s*$", f.read_text(), re.S).group(1))
    lig = d["ligand"]
    mol = Chem.MolFromMolBlock(lig["molblock"], sanitize=False)
    n, perms = mol.GetNumAtoms(), lig["automorphisms"]
    el = [a.GetSymbol() for a in mol.GetAtoms()]
    bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    ok = (list(range(n)) in perms
          and all(el[i] == el[p[i]] for p in perms for i in range(n))
          and all(tuple(sorted((p[a], p[b]))) in bonds for p in perms for a, b in bonds)
          and [a["name"] for a in lig["atoms"]] == [a["name"] for a in lig["atoms"]])
    failures += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {d['id']:<22} {len(perms):>3} automorphisms over {n} atoms")
raise SystemExit(failures)
