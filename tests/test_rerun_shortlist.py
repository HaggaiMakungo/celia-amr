import csv

from molecular_docking.rerun_shortlist import load_shortlist, write_shortlist


def test_shortlist_adds_and_persists_second_pass_columns(tmp_path):
    path = tmp_path / "shortlist.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("chembl_id", "pref_name"))
        writer.writeheader()
        writer.writerow({"chembl_id": "CHEMBL1", "pref_name": "Example"})
    rows, columns = load_shortlist(path)
    rows[0]["second_pass_status"] = "ok"
    write_shortlist(path, rows, columns)
    with path.open(newline="", encoding="utf-8") as handle:
        assert next(csv.DictReader(handle))["second_pass_status"] == "ok"
