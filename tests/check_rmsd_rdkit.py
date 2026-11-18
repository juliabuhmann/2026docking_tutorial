"""Cross-check the JS RMSDs (tests/out/*_random_poses.json, from test_metrics.js) against RDKit.

Ported from the docking_game prototype. Uses rdMolAlign.CalcRMS (no superposition) with our
bond-order-agnostic automorphisms passed as `map`. RDKit's own symmetry detection respects bond
orders and so can miss resonance-equivalent atoms (e.g. benzamidine's two amidinium N); that
value is printed for information only.
"""
import json
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
TOL_A = 1e-3


def load(name):
    text = (DATA / name).read_text()
    return json.loads(re.search(r"PW\.data\([^,]+,\s*[^,]+,\s*(.*)\);\s*$", text, re.S).group(1))


failures = 0
for dc in load("index.js")["docking"]["cases"]:
    cs = load(f"prot_{dc['level']}.js")
    ref = Chem.MolFromMolBlock(cs["ligand"]["molblock"], sanitize=False)
    poses = json.loads((ROOT / f"tests/out/{dc['level']}_random_poses.json").read_text())
    worst = rdkit_default_gap = 0.0
    maps = [[(i, j) for i, j in enumerate(perm)] for perm in cs["ligand"]["automorphisms"]]
    for p in poses:
        prb = Chem.Mol(ref)
        conf = prb.GetConformer()
        for i, xyz in enumerate(p["xyz"]):
            conf.SetAtomPosition(i, Point3D(*xyz))
        worst = max(worst, abs(rdMolAlign.CalcRMS(prb, ref, map=maps) - p["sym_rmsd"]))
        rdkit_default_gap = max(rdkit_default_gap, rdMolAlign.CalcRMS(prb, ref) - p["sym_rmsd"])
    ok = worst < TOL_A
    failures += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {dc['level']}: max |JS sym RMSD - RDKit CalcRMS| = {worst:.2e} A over {len(poses)} poses"
          f"  (RDKit default symmetry: up to {rdkit_default_gap:.2f} A higher)")
raise SystemExit(failures)
