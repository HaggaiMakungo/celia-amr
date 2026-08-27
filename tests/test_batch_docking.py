import csv

from molecular_docking.batch_docking import (
    DOCKING_COLUMNS,
    DockingOutcome,
    eligible_row_indices,
    load_manifest,
    run_batch,
    select_test_batch,
    write_manifest,
)


def test_eligible_rows_only_include_prepared_unattempted_ligands():
    rows = [
        {"prep_status": "ok", "docking_status": ""},
        {"prep_status": "ok", "docking_status": "not_attempted"},
        {"prep_status": "ok", "docking_status": "ok"},
        {"prep_status": "ok", "docking_status": "failed"},
        {"prep_status": "rejected_inorganic", "docking_status": ""},
    ]
    assert eligible_row_indices(rows) == [0, 1]


def test_test_batch_spans_the_molecular_weight_range():
    rows = [{"molecular_weight": str(weight)} for weight in range(100, 1100, 100)]
    selected = select_test_batch(list(range(len(rows))), rows, 4)
    assert selected[0] == 0
    assert selected[-1] == len(rows) - 1
    assert len(selected) == 4


def test_manifest_adds_and_persists_docking_columns(tmp_path):
    manifest = tmp_path / "ligand_library_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("chembl_id", "prep_status", "pdbqt_path"))
        writer.writeheader()
        writer.writerow({"chembl_id": "CHEMBL1", "prep_status": "ok", "pdbqt_path": "pdbqt/CHEMBL1.pdbqt"})
    rows, columns = load_manifest(manifest)
    assert set(DOCKING_COLUMNS).issubset(columns)
    rows[0]["docking_status"] = "ok"
    write_manifest(manifest, rows, columns)
    with manifest.open(newline="", encoding="utf-8") as handle:
        persisted = next(csv.DictReader(handle))
    assert persisted["docking_status"] == "ok"


def test_batch_updates_manifest_and_writes_failure_audit(monkeypatch, tmp_path):
    manifest = tmp_path / "ligand_library_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("chembl_id", "pref_name", "prep_status", "pdbqt_path", "molecular_weight"))
        writer.writeheader()
        writer.writerow({"chembl_id": "CHEMBL1", "pref_name": "Success", "prep_status": "ok", "pdbqt_path": "pdbqt/CHEMBL1.pdbqt", "molecular_weight": "100"})
        writer.writerow({"chembl_id": "CHEMBL2", "pref_name": "Failure", "prep_status": "ok", "pdbqt_path": "pdbqt/CHEMBL2.pdbqt", "molecular_weight": "500"})

    monkeypatch.setattr("molecular_docking.batch_docking.prepare_katg_target", lambda: object())

    def fake_dock(row_index, row, *_args):
        if row["chembl_id"] == "CHEMBL1":
            return DockingOutcome(row_index, "ok", "docking_results/CHEMBL1.pdbqt", "-7.000", attempted=True)
        return DockingOutcome(row_index, "failed", notes="simulated Vina failure", attempted=True)

    monkeypatch.setattr("molecular_docking.batch_docking._dock_one", fake_dock)
    summary = run_batch(tmp_path, test_batch_size=2, workers=1)
    rows, _ = load_manifest(manifest)
    assert summary == {"selected": 2, "ok": 1, "failed": 1, "not_attempted": 0}
    assert [row["docking_status"] for row in rows] == ["ok", "failed"]
    assert (tmp_path / "docking_failures.csv").read_text().count("simulated Vina failure") == 1
    assert "test batch (2 ligands)" in (tmp_path / "docking_summary.md").read_text()
