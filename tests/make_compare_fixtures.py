#!/usr/bin/env python
"""Fake "predictions" for testing docs/js/compare.js, written to tests/out/.

For each co-folding target the crystal chain and ligand are renumbered the way a structure
predictor numbers them (1..N along the deposited sequence), moved by a random rigid-body
transform, and written as
  <id>_exact.pdb   LocalFold-style PDB (protein chain A, ligand HETATM chain B)   expect ~0 A
  <id>_exact.cif   mmCIF, as AlphaFold Server writes                               expect ~0 A
  <id>_shift3.pdb  ligand moved 3 A before the transform                           expect 3.00 A
  <id>_renamed.pdb ligand atoms renamed L1, L2, ... (as after folding from SMILES)  expect a warning
Run with the `docking` env after scripts/build_data.py:  python tests/make_compare_fixtures.py
"""
import json
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data_cache"
OUT = ROOT / "tests" / "out"
SHIFT_A = 3.0
SEED = 7


def random_rotation(rng):
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def build(target, rng, shift=0.0, rename=False):
    st = gemmi.read_structure(str(CACHE / f"{target['pdb_id']}.cif"))
    st.setup_entities()
    st.remove_hydrogens()
    model = st[0]
    poly = model[target["chain"]].get_polymer()
    R, t = random_rotation(rng), rng.normal(scale=20, size=3)
    move = lambda p: R @ np.asarray(p) + t
    out = gemmi.Structure()
    m = gemmi.Model("1")
    prot = gemmi.Chain("A")
    for r in poly:
        if r.label_seq is None:
            continue
        nr = gemmi.Residue(); nr.name = r.name; nr.seqid = gemmi.SeqId(r.label_seq, " "); nr.entity_type = gemmi.EntityType.Polymer
        for a in r:
            if a.altloc not in ("\0", "A"):
                continue
            na = gemmi.Atom(); na.name = a.name; na.element = a.element; na.pos = gemmi.Position(*move(a.pos.tolist())); na.occ = 1; na.b_iso = 90
            nr.add_atom(na)
        prot.add_residue(nr)
    m.add_chain(prot)
    lig_chain = gemmi.Chain("B")
    lig = next(r for ch in model for r in ch if r.name == target["ligand_id"] and ch.name == target["chain"])
    nr = gemmi.Residue(); nr.name = lig.name; nr.seqid = gemmi.SeqId(1, " "); nr.het_flag = "H"; nr.entity_type = gemmi.EntityType.NonPolymer
    for k, a in enumerate(a for a in lig if a.altloc in ("\0", "A")):
        na = gemmi.Atom(); na.name = f"L{k + 1}" if rename else a.name; na.element = a.element
        na.pos = gemmi.Position(*move(np.array(a.pos.tolist()) + [shift, 0, 0])); na.occ = 1; na.b_iso = 90
        nr.add_atom(na)
    lig_chain.add_residue(nr)
    m.add_chain(lig_chain)
    out.add_model(m)
    out.setup_entities()
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    targets = json.loads((ROOT / "configs" / "cofold.json").read_text())["targets"]
    expect = {}
    for tg in targets:
        tid = tg["id"]
        build(tg, rng).write_pdb(str(OUT / f"{tid}_exact.pdb"))
        doc = build(tg, rng).make_mmcif_document()
        doc.write_file(str(OUT / f"{tid}_exact.cif"))
        build(tg, rng, shift=SHIFT_A).write_pdb(str(OUT / f"{tid}_shift3.pdb"))
        build(tg, rng, rename=True).write_pdb(str(OUT / f"{tid}_renamed.pdb"))
        expect[tid] = {"exact.pdb": 0.0, "exact.cif": 0.0, "shift3.pdb": SHIFT_A, "renamed.pdb": None}
    (OUT / "compare_expect.json").write_text(json.dumps(expect, indent=1))
    print(f"wrote fixtures for {len(targets)} targets to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
