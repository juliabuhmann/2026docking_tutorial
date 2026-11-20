#!/usr/bin/env python
"""Build every data file the Pocketworks site loads, from RCSB.

Inputs are the four files in configs/. Nothing structural is typed in by hand: PDB entries,
ligand definitions (CCD), sequences and citations are all downloaded here.

    configs/proteins.json  -> docs/data/prot_<id>.js   pocket-hunt levels (also used by the
                                                         algorithm and docking chapters)
    configs/morphs.json    -> docs/data/morph_<id>.js  apo/holo pairs, aligned and superposed
    configs/cofold.json    -> docs/data/cofold_<id>.js and docs/data/jobs/<id>_job.json
    configs/docking.json   -> copied into docs/data/index.js
                           -> docs/data/index.js       small catalogue the menu loads first

The data files are plain <script>s (`PW.data(...)`), not JSON, so the site also works when
index.html is opened straight from disk, where browsers block fetch().

Run with the workspace's `docking` conda env (rdkit, gemmi, scipy, numpy):
    python scripts/build_data.py
"""
import json
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"
CACHE = ROOT / "data_cache"
OUT = ROOT / "docs" / "data"

PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
CCD_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

WATER_NAMES = {"HOH", "WAT", "DOD"}
POLAR_ELEMENTS = {"N", "O"}
ACIDIC_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}

# Alpha spheres are shipped only inside this radius window, to keep the files small. The page's
# radius sliders stay inside it. Fpocket's own defaults (3.4 to 6.2 A) sit well inside.
ALPHA_SHIP_MIN_A = 2.8
ALPHA_SHIP_MAX_A = 8.0


# ----------------------------------------------------------------------------- downloads
def fetch(url: str, cache_name: str) -> str:
    CACHE.mkdir(exist_ok=True)
    path = CACHE / cache_name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as r:
            path.write_bytes(r.read())
    return path.read_text()


def fetch_json(url: str, cache_name: str) -> dict:
    return json.loads(fetch(url, cache_name))


def primary_citation(pdb_id: str) -> dict:
    entry = fetch_json(ENTRY_URL.format(pdb_id=pdb_id), f"{pdb_id}_entry.json")
    cit = next((c for c in entry.get("citation", []) if c.get("id") == "primary"), {})
    authors = cit.get("rcsb_authors") or []
    first = authors[0].split(",")[0] if authors else ""
    return {
        "title": entry["struct"]["title"],
        "resolution": (entry.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0],
        "citation": {
            "authors": f"{first} et al." if len(authors) > 1 else first,
            "journal": cit.get("rcsb_journal_abbrev") or cit.get("journal_abbrev") or "",
            "year": cit.get("year"),
            "title": cit.get("title") or "",
            "doi": cit.get("pdbx_database_id_doi"),
        },
    }


# ----------------------------------------------------------------------------- PDB parsing
def parse_atom_line(line: str) -> dict:
    return {
        "record": line[:6].strip(),
        "name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resn": line[17:20].strip(),
        "chain": line[21].strip(),
        "resi": int(line[22:26]),
        "icode": line[26].strip(),
        "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        "elem": (line[76:78].strip() or line[12:14].strip()).capitalize(),
        "line": line.rstrip("\n")[:80],
    }


def read_atoms(pdb_text: str):
    """All heavy atoms of the first model; alternate locations other than blank/A dropped."""
    atoms = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        a = parse_atom_line(line)
        if a["altloc"] not in ("", "A") or a["elem"] in ("H", "D"):
            continue
        atoms.append(a)
    return atoms


def first_instance(atoms, resn, chain=None):
    sel = [a for a in atoms if a["resn"] == resn and (chain is None or a["chain"] == chain)]
    if not sel:
        return []
    key = (sel[0]["chain"], sel[0]["resi"])
    return [a for a in sel if (a["chain"], a["resi"]) == key]


def pdb_block(atoms) -> str:
    # Altloc column blanked: the page should never see two copies of an atom.
    return "\n".join(a["line"][:16] + " " + a["line"][17:] for a in atoms) + "\nEND\n"


# ----------------------------------------------------------------------------- ligand chemistry
def build_ligand_mol(ligand_atoms, smiles_options):
    """Connectivity from the crystal geometry, bond orders from the CCD SMILES. The stereo SMILES
    is tried first; some (e.g. benzamidine's, with an explicit [H] on a C=N) do not match the
    hydrogen-free crystal ligand, so the plain SMILES is the fallback."""
    block = "\n".join(a["line"] for a in ligand_atoms) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False, proximityBonding=False)
    rdDetermineBonds.DetermineConnectivity(mol)
    for smi in smiles_options:
        try:
            return AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(Chem.MolFromSmiles(smi)), mol)
        except ValueError:
            continue
    print("    warning: bond orders not assigned; single bonds kept")
    return mol


def automorphisms(mol, cap=2000):
    """Heavy-atom graph automorphisms, bond orders ignored (so resonance-equivalent atoms such
    as the two amidinium N count as equivalent). Used for symmetry-corrected RMSD."""
    graph = Chem.RWMol(mol)
    for b in graph.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in graph.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
    # The query is written as SMARTS so that only element and connectivity are compared. Writing
    # SMARTS can renumber the atoms: query atom k is graph atom order[k]. Each match is mapped back
    # to graph numbering, perm[i] = the atom that atom i is matched to. (The prototype skipped this
    # step, which only went unnoticed because benzamidine's SMARTS keeps the original order.)
    smarts = Chem.MolToSmarts(graph)
    order = list(graph.GetProp("_smilesAtomOutputOrder").strip("[],").split(","))
    order = [int(x) for x in order if x != ""]
    query = Chem.MolFromSmarts(smarts)
    perms = set()
    for match in graph.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=cap):
        perm = [0] * graph.GetNumAtoms()
        for k, atom in enumerate(match):
            perm[order[k]] = atom
        perms.add(tuple(perm))
    perms = sorted(perms)
    assert tuple(range(graph.GetNumAtoms())) in perms, "identity missing from automorphisms"
    return perms


def polar_contacts(ligand_atoms, partners_all, cutoff):
    partners = [a for a in partners_all if a["elem"] in POLAR_ELEMENTS]
    pxyz = np.array([a["xyz"] for a in partners])
    out = []
    for i, la in enumerate(ligand_atoms):
        if la["elem"] not in POLAR_ELEMENTS:
            continue
        d = np.linalg.norm(pxyz - np.array(la["xyz"]), axis=1)
        for j in np.where(d <= cutoff)[0]:
            p = partners[j]
            kind = ("water" if p["resn"] in WATER_NAMES
                    else "salt bridge" if (p["resn"], p["name"]) in ACIDIC_ATOMS else "H-bond")
            out.append({"lig_index": i, "lig_atom": la["name"], "partner": f"{p['resn']}{p['resi']} {p['name']}",
                        "partner_resn": p["resn"], "partner_resi": p["resi"], "partner_chain": p["chain"],
                        "partner_xyz": p["xyz"], "distance": round(float(d[j]), 2), "kind": kind})
    return sorted(out, key=lambda c: (c["kind"] == "water", c["distance"]))


def residues_near(points, receptor, cutoff):
    pts = np.asarray(points, dtype=float)
    near = set()
    for a in receptor:
        if np.min(np.linalg.norm(pts - np.array(a["xyz"]), axis=1)) <= cutoff:
            near.add((a["chain"], a["resi"], a["resn"]))
    return [{"chain": c, "resi": r, "resn": n} for c, r, n in sorted(near)]


# ----------------------------------------------------------------------------- alpha spheres
def alpha_spheres(xyz: np.ndarray):
    """Circumspheres of the Delaunay tetrahedra of the atom centres (the alpha spheres of
    fpocket, Le Guilloux et al. 2009). Returns centres, radii and the 4 atom indices."""
    tri = Delaunay(xyz)
    simp = tri.simplices
    p = xyz[simp]                                  # (n, 4, 3)
    a = p[:, 1:] - p[:, :1]                        # (n, 3, 3)
    b = 0.5 * np.sum(a * a, axis=2)                # (n, 3)
    det = np.linalg.det(a)
    ok = np.abs(det) > 1e-6                        # skip flat tetrahedra
    c = np.full((len(simp), 3), np.nan)
    c[ok] = np.linalg.solve(a[ok], b[ok][..., None])[..., 0] + p[ok, 0]
    r = np.linalg.norm(c - p[:, 0], axis=1)
    keep = ok & (r >= ALPHA_SHIP_MIN_A) & (r <= ALPHA_SHIP_MAX_A)
    return c[keep], r[keep], simp[keep]


# ----------------------------------------------------------------------------- one level
def build_level(level: dict, cfg: dict, dock_cfg: dict) -> dict:
    pid = level["pdb_id"]
    print(f"[{level['id']}] {pid} / {level['ligand_id']}")
    atoms = read_atoms(fetch(PDB_URL.format(pdb_id=pid), f"{pid}.pdb"))
    keep_het = set(level.get("keep_hetero", []))
    chains = set(level["receptor_chains"])
    receptor = [a for a in atoms if a["chain"] in chains and
                (a["record"] == "ATOM" or a["resn"] in keep_het)]
    waters = [a for a in atoms if a["resn"] in WATER_NAMES]
    lig = first_instance(atoms, level["ligand_id"], level["ligand_chain"])
    if not lig:
        raise SystemExit(f"ligand {level['ligand_id']} not found in chain {level['ligand_chain']}")

    ccd = fetch_json(CCD_URL.format(ligand_id=level["ligand_id"]), f"{level['ligand_id']}_ccd.json")
    desc = ccd["rcsb_chem_comp_descriptor"]
    smiles = desc.get("SMILES_stereo") or desc["SMILES"]
    mol = build_ligand_mol(lig, [smiles, desc["SMILES"]])
    perms = automorphisms(mol)
    cutoff = cfg["pocket_residue_cutoff_A"]
    pocket = residues_near([a["xyz"] for a in lig], receptor, cutoff)
    pocket_keys = {(r["chain"], r["resi"]) for r in pocket}

    decoys = []
    for d in level.get("decoys", []):
        if "from_ligand" in d:
            het = first_instance(atoms, d["from_ligand"])
            res = residues_near([a["xyz"] for a in het], receptor, cutoff)
            centre = np.mean([a["xyz"] for a in het], axis=0)
            how = f"residues within {cutoff} A of {d['from_ligand']} in {pid}"
        else:
            res = [{"chain": d["chain"], "resi": r, "resn": next(
                (a["resn"] for a in receptor if a["chain"] == d["chain"] and a["resi"] == r), "?")}
                for r in d["residues"]]
            centre = np.mean([a["xyz"] for a in receptor
                              if a["chain"] == d["chain"] and a["resi"] in d["residues"]], axis=0)
            how = d["source"]
        overlap = [r for r in res if (r["chain"], r["resi"]) in pocket_keys]
        if overlap:
            print(f"    decoy '{d['name']}' shares {len(overlap)} residues with the true pocket; "
                  "those count as correct")
        res = [r for r in res if (r["chain"], r["resi"]) not in pocket_keys]
        decoys.append({"name": d["name"], "feedback": d["feedback"], "residues": res,
                       "centre": centre.round(2).tolist(), "defined_by": how})

    rec_xyz = np.array([a["xyz"] for a in receptor])
    centres, radii, verts = alpha_spheres(rec_xyz)
    contacts = polar_contacts(lig, receptor + waters, dock_cfg["polar_contact_cutoff_A"])
    contact_waters = {(c["partner_chain"], c["partner_resi"]) for c in contacts if c["kind"] == "water"}
    water_atoms = [w for w in waters if (w["chain"], w["resi"]) in contact_waters]
    meta = primary_citation(pid)

    print(f"    {len(receptor)} receptor atoms, ligand {len(lig)} heavy atoms, {len(perms)} automorphisms, "
          f"{len(pocket)} pocket residues, {len(radii)} alpha spheres shipped, {len(contacts)} polar contacts")

    return {
        "id": level["id"], "pdb_id": pid, "ligand_id": level["ligand_id"],
        "ligand_name": ccd["chem_comp"]["name"], "ligand_label": level["ligand_label"], "ligand_smiles": smiles,
        **meta,
        "receptor_pdb": pdb_block(receptor),
        "water_pdb": pdb_block(water_atoms),
        "ligand": {
            "atoms": [{"name": a["name"], "elem": a["elem"], "xyz": a["xyz"]} for a in lig],
            "molblock": Chem.MolToMolBlock(mol),
            "automorphisms": [list(p) for p in perms],
            "contacts": contacts,
        },
        "pocket_residues": pocket,
        "pocket_residue_cutoff_A": cutoff,
        "decoys": decoys,
        "alpha": {
            "x": centres[:, 0].round(2).tolist(), "y": centres[:, 1].round(2).tolist(),
            "z": centres[:, 2].round(2).tolist(), "r": radii.round(2).tolist(),
            "v": verts.astype(int).ravel().tolist(),
            "ship_min_A": ALPHA_SHIP_MIN_A, "ship_max_A": ALPHA_SHIP_MAX_A,
        },
    }


# ----------------------------------------------------------------------------- morphs
def modified_residues(pdb_text: str) -> dict:
    """MODRES records: (chain, resi, icode) -> (modified name, standard parent), e.g. PTR -> TYR."""
    out = {}
    for l in pdb_text.splitlines():
        if l.startswith("MODRES"):
            out[(l[16], int(l[18:22]), l[22].strip())] = (l[12:15].strip(), l[24:27].strip())
    return out


def chain_residues(atoms, chain, modres=None):
    """Residues of one chain. Modified residues listed in MODRES (HETATM records, e.g. a
    phosphotyrosine) are kept and named by their standard parent; 'modified' keeps the original."""
    res, modres = {}, modres or {}
    for a in atoms:
        if a["chain"] != chain:
            continue
        mod = modres.get((a["chain"], a["resi"], a["icode"]))
        if a["record"] != "ATOM" and not (mod and mod[0] == a["resn"]):
            continue
        r = res.setdefault((a["resi"], a["icode"]), {"resn": mod[1] if mod else a["resn"], "modified": mod[0] if mod else None, "atoms": {}})
        r["atoms"][a["name"]] = a
    return [v | {"key": k} for k, v in sorted(res.items())]


def one_letter(resn):
    info = gemmi.find_tabulated_residue(resn)
    return info.one_letter_code.upper() if info and info.one_letter_code != " " else "X"


def kabsch(P, Q):
    """R, t minimising |R P + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def build_morph(pair: dict, trim: float) -> dict:
    print(f"[morph {pair['id']}] {pair['start']['pdb_id']} -> {pair['end']['pdb_id']}")
    s_text = fetch(PDB_URL.format(pdb_id=pair["start"]["pdb_id"]), f"{pair['start']['pdb_id']}.pdb")
    e_text = fetch(PDB_URL.format(pdb_id=pair["end"]["pdb_id"]), f"{pair['end']['pdb_id']}.pdb")
    s_atoms, e_atoms = read_atoms(s_text), read_atoms(e_text)
    S = chain_residues(s_atoms, pair["start"]["chain"], modified_residues(s_text))
    E = chain_residues(e_atoms, pair["end"]["chain"], modified_residues(e_text))
    for tag, R in (("start", S), ("end", E)):
        for r in R:
            if r["modified"]:
                print(f"    {tag}: {r['modified']} {r['key'][0]} counted as {r['resn']}")
    s_seq = "".join(one_letter(r["resn"]) for r in S)
    e_seq = "".join(one_letter(r["resn"]) for r in E)
    aln = gemmi.align_string_sequences(list(s_seq), list(e_seq), [])
    pairs, i, j = [], 0, 0
    # walk the alignment strings ('-' marks a gap)
    a1, a2 = aln.add_gaps(s_seq, 1), aln.add_gaps(e_seq, 2)
    for c1, c2 in zip(a1, a2):
        if c1 != "-" and c2 != "-":
            pairs.append((i, j))
        i += c1 != "-"
        j += c2 != "-"
    print(f"    sequence identity {aln.calculate_identity():.1f}%, {len(pairs)} aligned residues")

    # common atoms: all names for identical residues, backbone + CB otherwise
    lines, s_xyz, e_xyz, ca_idx = [], [], [], []
    for i, j in pairs:
        rs, re_ = S[i], E[j]
        names = (list(re_["atoms"]) if rs["resn"] == re_["resn"]
                 else [n for n in ("N", "CA", "C", "O", "CB") if n in re_["atoms"]])
        for n in names:
            if n in rs["atoms"] and n in re_["atoms"]:
                if n == "CA":
                    ca_idx.append(len(lines))
                lines.append(re_["atoms"][n])
                s_xyz.append(rs["atoms"][n]["xyz"])
                e_xyz.append(re_["atoms"][n]["xyz"])
    s_xyz, e_xyz = np.array(s_xyz), np.array(e_xyz)

    # superpose start onto end on the CA atoms of the rigid core (iterative trimming)
    # Each round refits on the CA atoms within `trim` of the current fit. When a domain motion
    # leaves too few that close (adenylate kinase), keep the best-fitting 40% instead and refit.
    ca_all = np.array(ca_idx)
    core = ca_all
    for _ in range(30):
        R, t = kabsch(s_xyz[core], e_xyz[core])
        dev = np.linalg.norm(s_xyz[ca_all] @ R.T + t - e_xyz[ca_all], axis=1)
        new = ca_all[dev < trim]
        if len(new) < 0.4 * len(ca_all):
            new = ca_all[np.sort(np.argsort(dev)[: int(0.4 * len(ca_all))])]
        if np.array_equal(new, core):
            break
        core = new
    s_fit = s_xyz @ R.T + t
    ca_dev = np.linalg.norm(s_fit[ca_idx] - e_xyz[ca_idx], axis=1)
    print(f"    core {len(core)}/{len(ca_idx)} CA within {trim} A; overall CA RMSD "
          f"{np.sqrt((ca_dev ** 2).mean()):.2f} A, max CA shift {ca_dev.max():.1f} A")

    lig = []
    for code in pair["end"].get("ligand_ids", []):
        lig_chain = pair["end"].get("ligand_chain", pair["end"]["chain"])
        lig += [a for a in e_atoms if a["resn"] == code and a["chain"] == lig_chain]
    if pair["end"].get("ligand_ids") and not lig:
        raise SystemExit(f"ligand {pair['end']['ligand_ids']} not found in {pair['end']['pdb_id']}")
    print(f"    ligand: {len(lig)} heavy atoms")
    return {
        "id": pair["id"], "title": pair["title"], "story": pair["story"],
        "start": pair["start"] | primary_citation(pair["start"]["pdb_id"]),
        "end": pair["end"] | primary_citation(pair["end"]["pdb_id"]),
        "identity_pct": round(aln.calculate_identity(), 1),
        "core_ca": int(len(core)), "aligned_ca": len(ca_idx), "core_trim_A": trim,
        "ca_rmsd_A": round(float(np.sqrt((ca_dev ** 2).mean())), 2),
        "pdb": pdb_block(lines),            # atom names/residues of the end structure
        "start_xyz": s_fit.round(3).ravel().tolist(),
        "end_xyz": e_xyz.round(3).ravel().tolist(),
        "ligand_pdb": pdb_block(lig) if lig else "",
    }


# ----------------------------------------------------------------------------- co-folding
def build_cofold(target: dict, prot_data: dict) -> tuple[dict, list]:
    pid = target["pdb_id"]
    print(f"[cofold {target['id']}] {pid} chain {target['chain']} + {target['ligand_id']}")
    st = gemmi.read_structure(str(CACHE / f"{pid}.cif")) if (CACHE / f"{pid}.cif").exists() else None
    if st is None:
        fetch(CIF_URL.format(pdb_id=pid), f"{pid}.cif")
        st = gemmi.read_structure(str(CACHE / f"{pid}.cif"))
    st.setup_entities()
    ch = st[0][target["chain"]]
    poly = ch.get_polymer()
    ent = st.get_entity_of(poly)
    seq = gemmi.one_letter_code(ent.full_sequence).upper()
    # sequence position (1-based, as a structure predictor numbers it) -> deposited residue number
    seq_to_auth = {r.label_seq: r.seqid.num for r in poly if r.label_seq is not None}
    pocket_auth = {r["resi"] for r in prot_data["pocket_residues"] if r["chain"] == target["chain"]}
    ca = {}
    for r in poly:
        if r.label_seq is not None and r.find_atom("CA", "*"):
            ca[r.label_seq] = list(np.round(r.find_atom("CA", "*").pos.tolist(), 3))
    print(f"    {len(seq)} residues in SEQRES, {len(seq_to_auth)} modelled, "
          f"{sum(seq_to_auth[k] in pocket_auth for k in seq_to_auth)} pocket residues mapped")
    job = [{
        "name": f"{target['id']}_{pid}_{target['ligand_id']}",
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": 1}},
            {"ligand": {"ligand": f"CCD_{target['ligand_id']}", "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }]
    data = {
        "id": target["id"], "pdb_id": pid, "chain": target["chain"], "ligand_id": target["ligand_id"],
        "ligand_name": prot_data["ligand_name"], "ligand_label": prot_data["ligand_label"], "ligand_smiles": prot_data["ligand_smiles"],
        "note": target.get("note", ""), "sequence": seq,
        "ca_by_seq": {str(k): v for k, v in sorted(ca.items())},
        "pocket_seq": sorted(k for k, v in seq_to_auth.items() if v in pocket_auth),
        "seq_to_auth": {str(k): v for k, v in sorted(seq_to_auth.items())},
        "job_file": f"data/jobs/{target['id']}_job.json",
    }
    return data, job


# ----------------------------------------------------------------------------- main
def write_js(path: Path, kind: str, key: str, obj) -> None:
    path.write_text("// Generated by scripts/build_data.py from RCSB data. Do not edit by hand.\n"
                    f"PW.data({json.dumps(kind)}, {json.dumps(key)}, {json.dumps(obj, separators=(',', ':'))});\n")
    print(f"    wrote {path.relative_to(ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


def main():
    prot_cfg = json.loads((CONFIGS / "proteins.json").read_text())
    morph_cfg = json.loads((CONFIGS / "morphs.json").read_text())
    cofold_cfg = json.loads((CONFIGS / "cofold.json").read_text())
    dock_cfg = json.loads((CONFIGS / "docking.json").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "jobs").mkdir(exist_ok=True)

    index = {"levels": [], "morphs": [], "cofold": [], "docking": dock_cfg}
    prot = {}
    for level in prot_cfg["levels"]:
        d = build_level(level, prot_cfg, dock_cfg)
        prot[level["id"]] = d
        write_js(OUT / f"prot_{level['id']}.js", "prot", level["id"], d)
        index["levels"].append({k: level[k] for k in ("id", "pdb_id", "title", "difficulty", "story")} | {
            "caveat": level.get("caveat", ""), "story_refs": level.get("story_refs", ""), "intro": level.get("intro", []),
            "ligand_id": level["ligand_id"], "ligand_name": d["ligand_name"], "ligand_label": level["ligand_label"],
            "n_ligand_atoms": len(d["ligand"]["atoms"]), "chains": level["receptor_chains"],
            "resolution": d["resolution"], "citation": d["citation"]})

    for pair in morph_cfg["pairs"]:
        d = build_morph(pair, morph_cfg["core_trim_A"])
        write_js(OUT / f"morph_{pair['id']}.js", "morph", pair["id"], d)
        index["morphs"].append({k: d[k] for k in ("id", "title", "ca_rmsd_A")})

    for target in cofold_cfg["targets"]:
        d, job = build_cofold(target, prot[target["id"]])
        write_js(OUT / f"cofold_{target['id']}.js", "cofold", target["id"], d)
        (OUT / "jobs" / f"{target['id']}_job.json").write_text(json.dumps(job, indent=2) + "\n")
        index["cofold"].append({k: d[k] for k in ("id", "pdb_id", "ligand_id", "ligand_name", "ligand_label", "note")}
                               | {"length": len(d["sequence"])})

    (OUT / "index.js").write_text("// Generated by scripts/build_data.py. Do not edit by hand.\n"
                                  f"PW.data('index', 'index', {json.dumps(index, indent=1)});\n")
    print("wrote docs/data/index.js")


if __name__ == "__main__":
    main()
