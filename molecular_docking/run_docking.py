"""Branch 3: Molecular Docking.

Runs a real AutoDock Vina docking of a drug compound against a TB drug
target, for the same katG/Isoniazid and rpoB/Rifampicin pairs used by
`genomic_analysis`. No mocked scores: this fetches a real PDB structure and
a real PubChem compound, prepares both for docking, and calls the actual
Vina binary as a subprocess.

Standalone script — no dependency on the other three branches. Run directly:

    python run_docking.py --gene katG

## Why this doesn't use conda / the `vina` PyPI package / OpenBabel

`environment.yml` (and the PyPI `vina` package) assume a working conda/Boost
toolchain, which this machine doesn't have, and neither the `vina` nor
`openbabel`/`openbabel-wheel` PyPI packages ship Windows wheels — `vina`'s
sdist requires Boost to build from source. This module instead uses:

- The official precompiled `vina.exe` from the AutoDock-Vina GitHub releases
  (downloaded once into `molecular_docking/bin/`, invoked as a subprocess) —
  the same binary conda would have installed, just fetched directly.
- `meeko` (the AutoDock team's modern, pure-Python ligand/receptor PDBQT
  preparation tool, built on RDKit) in place of OpenBabel/AutoDockTools for
  file conversion.

If conda *is* available in your environment, `environment.yml`'s
`autodock-vina`/`openbabel` conda packages work as an equivalent substitute
for the binary download + meeko path used here.

## Docking targets

- **katG / Isoniazid**: PDB `2CCA` (M. tuberculosis KatG, X-ray, 2.0 Angstrom;
  notably this is the same structure that also contains the S315T resistance
  mutant relevant to `genomic_analysis`). KatG activates isoniazid at its
  heme active site, so the docking box is centered on the real heme (HEM)
  cofactor's coordinates, read directly from the downloaded structure.
  KatG carries a well-documented autocatalytic Met255-Tyr229-Trp107 (MYW)
  covalent crosslink (Yu et al. 2003; Yamada et al. 2001) essential for its
  catalase activity. meeko's default residue perception doesn't have a
  template for this crosslink and mis-parses it as a malformed inter-residue
  bond, so that specific bond is deleted before receptor parameterization
  (`KATG_BONDS_TO_DELETE`) — the crosslink itself is real and left in the
  coordinates; only its interpretation as a "residue-connecting" bond by the
  parser is removed, so the crosslinked residues dock as (very slightly)
  simplified rigid rather than a rejected structure.
- **rpoB / Rifampicin**: PDB `5UH6` (M. tuberculosis RNA polymerase -
  rifampicin cryo-EM complex). This is a large multi-subunit complex; it is
  wired into this module's target registry but has not been run end-to-end
  here (see README) — expect longer receptor-prep and docking runtimes.
"""

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import meeko
import truststore

truststore.inject_into_ssl()

import pubchempy as pcp  # noqa: E402  (must follow truststore.inject_into_ssl())
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parent
STRUCTURES_DIR = REPO_ROOT / "data" / "structures"
BIN_DIR = MODULE_DIR / "bin"
WORK_DIR = MODULE_DIR / "work"

VINA_EXE_PATH = BIN_DIR / "vina.exe"
VINA_DOWNLOAD_URL = (
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"
)


@dataclass(frozen=True)
class DockingTarget:
    gene: str
    drug_name: str
    pdb_id: str
    box_ligand_resname: str  # HETATM residue used to derive the docking box center
    box_chain: str
    bonds_to_delete: tuple  # [(res_id_a, res_id_b), ...] -- see module docstring
    strip_resnames: tuple  # HETATM residues removed before receptor prep

    @property
    def receptor_pdb_path(self) -> Path:
        return STRUCTURES_DIR / f"{self.pdb_id}.pdb"


TARGETS = {
    "katG": DockingTarget(
        gene="katG",
        drug_name="Isoniazid",
        pdb_id="2CCA",
        box_ligand_resname="HEM",
        box_chain="A",
        bonds_to_delete=(("A:229", "A:255"), ("B:229", "B:255")),
        strip_resnames=("HOH", "HEM"),
    ),
    "rpoB": DockingTarget(
        gene="rpoB",
        drug_name="Rifampicin",
        pdb_id="5UH6",
        box_ligand_resname="RFP",
        box_chain="C",
        bonds_to_delete=(),
        strip_resnames=("HOH",),
    ),
}


def download_pdb_structure(pdb_id: str, out_path: Path) -> None:
    import urllib.request

    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    with urllib.request.urlopen(url, timeout=60) as response:
        out_path.write_bytes(response.read())


def ensure_vina_binary(exe_path: Path = VINA_EXE_PATH) -> Path:
    if exe_path.exists():
        return exe_path
    import urllib.request

    exe_path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(VINA_DOWNLOAD_URL, timeout=120) as response:
        exe_path.write_bytes(response.read())
    return exe_path


def fetch_ligand_smiles(drug_name: str) -> tuple[int, str]:
    """Look up a compound on PubChem by name; returns (CID, isomeric SMILES)."""
    compounds = pcp.get_compounds(drug_name, "name")
    if not compounds:
        raise ValueError(f"No PubChem compound found for {drug_name!r}")
    compound = compounds[0]
    return compound.cid, compound.smiles


def build_ligand_sdf(smiles: str, name: str, out_path: Path, seed: int = 42) -> None:
    """Embed a real 3D conformer for a SMILES string and write it to SDF."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=seed)
    if embed_result != 0:
        raise RuntimeError(f"RDKit 3D embedding failed for {name}")
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", name)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()


def prepare_ligand_pdbqt(sdf_path: Path, out_pdbqt: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", str(sdf_path), "-o", str(out_pdbqt)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"meeko ligand prep failed:\n{result.stdout}\n{result.stderr}")


def clean_receptor_pdb(pdb_path: Path, strip_resnames: tuple, out_path: Path) -> None:
    lines = pdb_path.read_text().splitlines(keepends=True)
    kept = [
        line
        for line in lines
        if not (line.startswith("HETATM") and line[17:20].strip() in strip_resnames)
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(kept))


def prepare_receptor_pdbqt(cleaned_pdb_path: Path, bonds_to_delete: tuple, out_pdbqt: Path) -> None:
    """Parameterize a receptor and write its rigid PDBQT using meeko's Polymer API directly.

    Uses the Polymer API (rather than the `mk_prepare_receptor` CLI) because
    only the Python API exposes `bonds_to_delete`, needed to handle katG's
    MYW covalent crosslink (see module docstring).
    """
    chem_templates = meeko.ResidueChemTemplates.create_from_defaults()
    mk_prep = meeko.MoleculePreparation(load_atom_params="ad4_types", charge_model="gasteiger")
    pdb_string = cleaned_pdb_path.read_text()

    polymer = meeko.Polymer.from_pdb_string(
        pdb_string,
        chem_templates,
        mk_prep,
        allow_bad_res=True,
        bonds_to_delete=list(bonds_to_delete) if bonds_to_delete else None,
        default_altloc="A",
    )
    rigid_pdbqt, flex_dict = meeko.PDBQTWriterLegacy.write_from_polymer(polymer)

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.write_text(rigid_pdbqt)


def compute_box_center(pdb_path: Path, resname: str, chain: str) -> tuple[float, float, float]:
    """Real docking-box center: the centroid of a real cofactor/co-crystallized ligand."""
    coords = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == resname and line[21] == chain:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ValueError(f"No {resname!r} atoms found in chain {chain!r} of {pdb_path}")
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


@dataclass(frozen=True)
class DockingPose:
    mode: int
    affinity_kcal_per_mol: float
    rmsd_lower_bound: float
    rmsd_upper_bound: float


VINA_POSE_LINE = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.?\d*|0)\s+(\d+\.?\d*|0)\s*$"
)


def parse_vina_output(stdout: str) -> list[DockingPose]:
    poses = []
    for line in stdout.splitlines():
        match = VINA_POSE_LINE.match(line)
        if match:
            mode, affinity, rmsd_lb, rmsd_ub = match.groups()
            poses.append(
                DockingPose(
                    mode=int(mode),
                    affinity_kcal_per_mol=float(affinity),
                    rmsd_lower_bound=float(rmsd_lb),
                    rmsd_upper_bound=float(rmsd_ub),
                )
            )
    return poses


def run_vina(
    vina_exe: Path,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box_center: tuple,
    box_size: tuple,
    out_pdbqt: Path,
    seed: int = 42,
    exhaustiveness: int = 8,
    cpu: int | None = None,
    timeout_seconds: float | None = None,
) -> list[DockingPose]:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(vina_exe),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box_center[0]),
        "--center_y", str(box_center[1]),
        "--center_z", str(box_center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--out", str(out_pdbqt),
        "--seed", str(seed),
        "--exhaustiveness", str(exhaustiveness),
    ]
    if cpu is not None:
        cmd.extend(["--cpu", str(cpu)])
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"vina timed out after {timeout_seconds:.0f} seconds") from exc
    if result.returncode != 0:
        raise RuntimeError(f"vina failed:\n{result.stdout}\n{result.stderr}")
    return parse_vina_output(result.stdout)


def dock(gene: str, box_padding: float = 20.0, exhaustiveness: int = 8, seed: int = 42) -> dict:
    target = TARGETS[gene]
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    if not target.receptor_pdb_path.exists():
        download_pdb_structure(target.pdb_id, target.receptor_pdb_path)

    vina_exe = ensure_vina_binary()

    box_center = compute_box_center(target.receptor_pdb_path, target.box_ligand_resname, target.box_chain)

    cleaned_pdb = WORK_DIR / f"{target.pdb_id}_cleaned.pdb"
    clean_receptor_pdb(target.receptor_pdb_path, target.strip_resnames, cleaned_pdb)

    receptor_pdbqt = WORK_DIR / f"{target.pdb_id}_receptor.pdbqt"
    prepare_receptor_pdbqt(cleaned_pdb, target.bonds_to_delete, receptor_pdbqt)

    cid, smiles = fetch_ligand_smiles(target.drug_name)
    ligand_name = f"{target.drug_name.lower()}_CID{cid}"
    ligand_sdf = WORK_DIR / f"{ligand_name}.sdf"
    build_ligand_sdf(smiles, ligand_name, ligand_sdf, seed=seed)

    ligand_pdbqt = WORK_DIR / f"{ligand_name}.pdbqt"
    prepare_ligand_pdbqt(ligand_sdf, ligand_pdbqt)

    docked_out = WORK_DIR / f"{ligand_name}_docked.pdbqt"
    poses = run_vina(
        vina_exe,
        receptor_pdbqt,
        ligand_pdbqt,
        box_center,
        (box_padding, box_padding, box_padding),
        docked_out,
        seed=seed,
        exhaustiveness=exhaustiveness,
    )

    return {
        "gene": target.gene,
        "drug": target.drug_name,
        "pubchem_cid": cid,
        "pdb_id": target.pdb_id,
        "box_center": box_center,
        "box_size": (box_padding, box_padding, box_padding),
        "poses": poses,
        "best_affinity_kcal_per_mol": poses[0].affinity_kcal_per_mol if poses else None,
        "docked_pdbqt_path": docked_out,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="TB molecular docking (real AutoDock Vina run)")
    parser.add_argument("--gene", required=True, choices=sorted(TARGETS.keys()))
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--box-padding", type=float, default=20.0)
    args = parser.parse_args()

    result = dock(args.gene, box_padding=args.box_padding, exhaustiveness=args.exhaustiveness, seed=args.seed)

    print(f"Gene: {result['gene']}  Drug: {result['drug']} (PubChem CID {result['pubchem_cid']})")
    print(f"Receptor: PDB {result['pdb_id']}")
    cx, cy, cz = result["box_center"]
    print(f"Docking box center: ({cx:.3f}, {cy:.3f}, {cz:.3f})  size: {result['box_size']} A")
    print(f"Best affinity: {result['best_affinity_kcal_per_mol']:.3f} kcal/mol")
    print()
    print(" mode |   affinity | rmsd l.b. | rmsd u.b.")
    for pose in result["poses"]:
        print(f"  {pose.mode:>3} | {pose.affinity_kcal_per_mol:>10.3f} | {pose.rmsd_lower_bound:>9.3f} | {pose.rmsd_upper_bound:>9.3f}")
    print(f"\nDocked poses written to {result['docked_pdbqt_path']}")


if __name__ == "__main__":
    main()
