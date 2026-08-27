"""Build a versioned, docking-ready ChEMBL approved-drug ligand library.

This module deliberately reuses ``run_docking.build_ligand_sdf`` and
``run_docking.prepare_ligand_pdbqt``.  Those are the RDKit embedding,
MMFF optimisation, and Meeko PDBQT-conversion settings used for the verified
isoniazid/KatG baseline, so every batch ligand is prepared equivalently.

Run from the repository root with the supported Python 3.11/conda environment:

    python molecular_docking/build_ligand_library.py

The generated data are intentionally ignored by Git.  Re-run this command to
refresh the library against a newer ChEMBL release.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urljoin

import truststore
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.SaltRemover import SaltRemover

from molecular_docking.run_docking import build_ligand_sdf, prepare_ligand_pdbqt

truststore.inject_into_ssl()

try:
    import requests
except ImportError as exc:  # pragma: no cover - exercised only in a broken setup
    raise RuntimeError("This tool needs the project's requests dependency.") from exc


MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parent
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "ligand_library" / "chembl_approved_v1"
CHEMBL_API_ROOT = "https://www.ebi.ac.uk/chembl/api/data"
MANIFEST_COLUMNS = (
    "chembl_id",
    "pref_name",
    "canonical_smiles",
    "prepared_smiles",
    "molecular_weight",
    "pdbqt_path",
    "drugbank_id",
    "pubchem_cid",
    "prep_status",
    "prep_notes",
)


@dataclass(frozen=True)
class LigandRecord:
    chembl_id: str
    pref_name: str
    canonical_smiles: str
    molecule_type: str
    max_phase: int
    parent_chembl_id: str
    drugbank_id: str = ""
    pubchem_cid: str = ""
    prepared_smiles: str = ""
    molecular_weight: float | None = None
    pdbqt_path: str = ""
    prep_status: str = "pending"
    prep_notes: str = ""


def _response_json(response: Any) -> dict[str, Any]:
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise RuntimeError("ChEMBL returned an unexpected non-object JSON response")
    return payload


def fetch_approved_molecules(session: requests.Session, api_root: str = CHEMBL_API_ROOT) -> list[dict[str, Any]]:
    """Retrieve every ChEMBL phase-4 small-molecule record, following pages."""
    url = f"{api_root}/molecule.json"
    params: dict[str, Any] | None = {
        "max_phase": 4,
        "molecule_type": "Small molecule",
        "limit": 1000,
    }
    records: list[dict[str, Any]] = []
    while url:
        response = session.get(url, params=params, timeout=90)
        payload = _response_json(response)
        page_records = payload.get("molecules", [])
        if not isinstance(page_records, list):
            raise RuntimeError("ChEMBL response did not contain a molecule list")
        records.extend(record for record in page_records if isinstance(record, dict))
        next_url = payload.get("page_meta", {}).get("next")
        # ChEMBL currently returns paths such as
        # /chembl/api/data/molecule.json?... rather than absolute URLs.
        url = urljoin(api_root + "/", next_url) if isinstance(next_url, str) and next_url else ""
        params = None  # ChEMBL's next URL already has its cursor and filters.
    return records


def fetch_chembl_release(session: requests.Session, api_root: str = CHEMBL_API_ROOT) -> str:
    """Return the database release advertised by the ChEMBL status endpoint."""
    payload = _response_json(session.get(f"{api_root}/status.json", timeout=30))
    for key in ("chembl_db_version", "chembl_release", "version"):
        value = payload.get(key)
        if value:
            return str(value)
    return "unknown (status endpoint did not expose a release field)"


def _xref_ids(raw_record: dict[str, Any]) -> tuple[str, str]:
    """Extract optional DrugBank and PubChem IDs without assuming one API shape."""
    drugbank_id = ""
    pubchem_cid = ""
    for xref in raw_record.get("cross_references") or []:
        if not isinstance(xref, dict):
            continue
        source = " ".join(
            str(xref.get(key, "")) for key in ("xref_src", "source", "xref_type")
        ).lower()
        identifier = next(
            (str(xref[key]) for key in ("xref_id", "xref_name", "id") if xref.get(key) is not None),
            "",
        )
        if "drugbank" in source and not drugbank_id:
            drugbank_id = identifier
        if "pubchem" in source and not pubchem_cid:
            pubchem_cid = identifier.removeprefix("CID:").strip()
    return drugbank_id, pubchem_cid


def record_from_chembl(raw_record: dict[str, Any]) -> LigandRecord | None:
    """Normalize the fields used by the library, ignoring malformed API records."""
    chembl_id = str(raw_record.get("molecule_chembl_id") or "").strip()
    structures = raw_record.get("molecule_structures") or {}
    canonical_smiles = structures.get("canonical_smiles") if isinstance(structures, dict) else None
    if not chembl_id or not canonical_smiles:
        return None
    drugbank_id, pubchem_cid = _xref_ids(raw_record)
    hierarchy = raw_record.get("molecule_hierarchy") or {}
    parent = hierarchy.get("parent_chembl_id") if isinstance(hierarchy, dict) else None
    return LigandRecord(
        chembl_id=chembl_id,
        pref_name=str(raw_record.get("pref_name") or ""),
        canonical_smiles=str(canonical_smiles),
        molecule_type=str(raw_record.get("molecule_type") or ""),
        # ChEMBL 37 serializes this numeric field as the string "4.0".
        max_phase=int(float(raw_record.get("max_phase") or 0)),
        parent_chembl_id=str(parent or chembl_id),
        drugbank_id=drugbank_id,
        pubchem_cid=pubchem_cid,
    )


def _normalise_structure(record: LigandRecord, salt_remover: SaltRemover) -> LigandRecord:
    mol = Chem.MolFromSmiles(record.canonical_smiles)
    if mol is None:
        return replace(record, prep_status="rejected_unparseable", prep_notes="RDKit MolFromSmiles returned None")
    parent = salt_remover.StripMol(mol, dontRemoveEverything=True)
    prepared_smiles = Chem.MolToSmiles(parent, canonical=True)
    molecular_weight = Descriptors.MolWt(parent)
    return replace(record, prepared_smiles=prepared_smiles, molecular_weight=molecular_weight)


def _fragment_count(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return len(Chem.GetMolFrags(mol)) if mol is not None else 0


def _contains_carbon(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None and any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())


def _resolve_to_chembl_parent(record: LigandRecord, normalized_by_id: dict[str, LigandRecord]) -> LigandRecord:
    """Use ChEMBL's own parent relationship when RDKit misses an organic salt.

    RDKit's default SaltRemover does not recognise every organic counterion
    (besylate, xinafoate, pamoate, citrate, ...). ChEMBL's molecule hierarchy
    identifies the canonical active parent without guessing which fragment is
    pharmacologically relevant.
    """
    if (
        record.prep_status != "pending"
        or _fragment_count(record.prepared_smiles) <= 1
        or record.parent_chembl_id == record.chembl_id
    ):
        return record
    parent = normalized_by_id.get(record.parent_chembl_id)
    if parent is None or parent.prep_status != "pending":
        return replace(
            record,
            prep_status="rejected_unresolved_multifragment",
            prep_notes=(
                f"Multi-fragment structure; ChEMBL parent {record.parent_chembl_id} "
                "was not present as a usable record in this source pull"
            ),
        )
    return replace(
        record,
        prepared_smiles=parent.prepared_smiles,
        molecular_weight=parent.molecular_weight,
        prep_notes=f"Resolved multi-fragment salt/form to ChEMBL parent {parent.chembl_id}",
    )


def normalise_and_deduplicate(records: Iterable[LigandRecord], max_molecular_weight: float = 1500.0) -> list[LigandRecord]:
    """Normalise, clean, and retain one candidate per docking-ready structure.

    All source records remain in the returned list: non-selected salt/duplicate
    records are marked ``excluded_duplicate`` for a complete audit trail.
    """
    salt_remover = SaltRemover()
    initially_normalized = [_normalise_structure(record, salt_remover) for record in records]
    normalized_by_id = {record.chembl_id: record for record in initially_normalized}
    processed: list[LigandRecord] = []
    for record in initially_normalized:
        resolved = _resolve_to_chembl_parent(record, normalized_by_id)
        if resolved.prep_status == "pending" and not _contains_carbon(resolved.prepared_smiles):
            resolved = replace(
                resolved,
                prep_status="rejected_inorganic",
                prep_notes="Prepared structure contains no carbon atoms",
            )
        if resolved.prep_status == "pending" and _fragment_count(resolved.prepared_smiles) > 1:
            resolved = replace(
                resolved,
                prep_status="rejected_unresolved_multifragment",
                prep_notes=(
                    "Multi-fragment structure has no distinct usable ChEMBL parent; "
                    "manual review is required before docking"
                ),
            )
        processed.append(resolved)
    groups: dict[str, list[LigandRecord]] = {}
    for record in processed:
        if record.prep_status != "pending":
            continue
        # Salt removal makes the prepared canonical SMILES the reliable
        # structure-level key. It catches both ChEMBL parent/child entries and
        # separate hierarchy records that resolve to the same parent structure.
        group_key = record.prepared_smiles
        groups.setdefault(group_key, []).append(record)

    result: dict[str, LigandRecord] = {record.chembl_id: record for record in processed}
    for candidates in groups.values():
        selected = min(
            candidates,
            key=lambda record: (
                record.chembl_id != record.parent_chembl_id,
                record.molecular_weight if record.molecular_weight is not None else float("inf"),
                record.chembl_id,
            ),
        )
        if selected.molecular_weight is not None and selected.molecular_weight > max_molecular_weight:
            selected = replace(
                selected,
                prep_status="rejected_size",
                prep_notes=f"Prepared parent molecular weight {selected.molecular_weight:.2f} Da exceeds {max_molecular_weight:.0f} Da",
            )
        result[selected.chembl_id] = selected
        for record in candidates:
            if record.chembl_id != selected.chembl_id:
                result[record.chembl_id] = replace(
                    record,
                    prep_status="excluded_duplicate",
                    prep_notes=f"Duplicate/salt form; selected {selected.chembl_id} for this parent structure",
                )
    return sorted(result.values(), key=lambda record: record.chembl_id)


def _valid_pdbqt(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0 and "ROOT" in path.read_text(errors="replace")


def prepare_library(records: Iterable[LigandRecord], output_dir: Path, seed: int = 42) -> list[LigandRecord]:
    """Make one PDBQT for every eligible record using the isoniazid baseline path."""
    pdbqt_dir = output_dir / "pdbqt"
    staging_dir = output_dir / ".staging_sdf"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    staging_dir.mkdir(parents=True, exist_ok=True)
    prepared: list[LigandRecord] = []
    try:
        for record in records:
            if record.prep_status != "pending":
                prepared.append(record)
                continue
            sdf_path = staging_dir / f"{record.chembl_id}.sdf"
            pdbqt_path = pdbqt_dir / f"{record.chembl_id}.pdbqt"
            try:
                build_ligand_sdf(record.prepared_smiles, record.chembl_id, sdf_path, seed=seed)
                prepare_ligand_pdbqt(sdf_path, pdbqt_path)
                if not _valid_pdbqt(pdbqt_path):
                    raise RuntimeError("Meeko created an empty or invalid PDBQT")
                prepared.append(
                    replace(record, pdbqt_path=pdbqt_path.relative_to(output_dir).as_posix(), prep_status="ok")
                )
            except Exception as exc:  # keep processing and make failures auditable
                if pdbqt_path.exists():
                    pdbqt_path.unlink()
                prepared.append(replace(record, prep_status="failed_preparation", prep_notes=str(exc)))
            finally:
                if sdf_path.exists():
                    sdf_path.unlink()
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)
    return prepared


def write_manifest(records: Iterable[LigandRecord], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "ligand_library_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "chembl_id": record.chembl_id,
                    "pref_name": record.pref_name,
                    "canonical_smiles": record.canonical_smiles,
                    "prepared_smiles": record.prepared_smiles,
                    "molecular_weight": "" if record.molecular_weight is None else f"{record.molecular_weight:.4f}",
                    "pdbqt_path": record.pdbqt_path,
                    "drugbank_id": record.drugbank_id,
                    "pubchem_cid": record.pubchem_cid,
                    "prep_status": record.prep_status,
                    "prep_notes": record.prep_notes,
                }
            )
    return manifest_path


def write_source_metadata(output_dir: Path, release: str, raw_count: int, retained_count: int) -> Path:
    source_path = output_dir / "SOURCE.md"
    source_path.write_text(
        "# ChEMBL approved-drug ligand library source\n\n"
        f"- ChEMBL database release: `{release}`\n"
        f"- Retrieved (UTC): `{datetime.now(UTC).isoformat()}`\n"
        "- Endpoint: `https://www.ebi.ac.uk/chembl/api/data/molecule.json`\n"
        "- Query filters: `max_phase=4`, `molecule_type=Small molecule`\n"
        f"- Raw queried records: `{raw_count}`\n"
        f"- Records with a canonical SMILES: `{retained_count}`\n"
        "- Preparation: RDKit ETKDG embedding (seed 42), MMFF optimisation, and "
        "Meeko `mk_prepare_ligand`; these are imported directly from "
        "`molecular_docking/run_docking.py`, the verified isoniazid baseline.\n"
        "- Protonation: no separate pH adjustment is applied because the baseline "
        "pipeline did not apply one; this preserves score comparability.\n"
        "- Cleanup: RDKit salt removal is applied to every record. Remaining "
        "multi-fragment salts are resolved to an already-pulled ChEMBL parent where "
        "available; otherwise they are flagged for review. Carbon-free prepared "
        "structures are rejected as inorganic.\n"
        "- Deduplication: select one smallest representative per canonical parent "
        "structure after cleanup; excluded source records remain in the manifest.\n",
        encoding="utf-8",
    )
    return source_path


def build_library(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    expected_minimum: int = 1500,
    expected_maximum: int = 2500,
    allow_unexpected_record_count: bool = False,
    seed: int = 42,
) -> tuple[list[LigandRecord], str]:
    """Fetch, validate, prepare, and document a complete ChEMBL library."""
    with requests.Session() as session:
        raw_records = fetch_approved_molecules(session)
        release = fetch_chembl_release(session)
    raw_count = len(raw_records)
    records = [record for raw in raw_records if (record := record_from_chembl(raw)) is not None]
    filtered = normalise_and_deduplicate(records)
    ready_count = sum(record.prep_status == "pending" for record in filtered)
    if not allow_unexpected_record_count and not expected_minimum <= ready_count <= expected_maximum:
        raise RuntimeError(
            f"Post-cleanup ChEMBL library contains {ready_count} docking-ready records, "
            f"outside the expected {expected_minimum}-{expected_maximum} range. Stopping "
            "before preparation; inspect the cleanup result or rerun with "
            "--allow-unexpected-record-count after review."
        )
    completed = prepare_library(filtered, output_dir, seed=seed)
    write_manifest(completed, output_dir)
    write_source_metadata(output_dir, release, raw_count, len(records))
    return completed, release


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the ChEMBL approved-drug PDBQT library")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--expected-minimum", type=int, default=1500)
    parser.add_argument("--expected-maximum", type=int, default=2500)
    parser.add_argument("--allow-unexpected-record-count", action="store_true")
    args = parser.parse_args()
    records, release = build_library(
        output_dir=args.output_dir,
        expected_minimum=args.expected_minimum,
        expected_maximum=args.expected_maximum,
        allow_unexpected_record_count=args.allow_unexpected_record_count,
        seed=args.seed,
    )
    counts: dict[str, int] = {}
    for record in records:
        counts[record.prep_status] = counts.get(record.prep_status, 0) + 1
    print(f"ChEMBL release: {release}")
    print(f"Manifest: {args.output_dir / 'ligand_library_manifest.csv'}")
    print("Status counts: " + json.dumps(counts, sort_keys=True))
    if counts.get("ok", 0) == 0:
        sys.exit("No ligands prepared successfully; inspect the manifest before docking.")


if __name__ == "__main__":
    main()
