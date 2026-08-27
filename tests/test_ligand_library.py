from pathlib import Path

from molecular_docking.build_ligand_library import (
    LigandRecord,
    fetch_approved_molecules,
    normalise_and_deduplicate,
    prepare_library,
    record_from_chembl,
    write_manifest,
)


def test_fetch_approved_molecules_follows_relative_next_page_urls():
    class Response:
        def __init__(self, payload):
            self.payload = payload

        def raise_for_status(self):
            return None

        def json(self):
            return self.payload

    class Session:
        def __init__(self):
            self.urls = []

        def get(self, url, params=None, timeout=None):
            self.urls.append((url, params))
            if len(self.urls) == 1:
                return Response({"molecules": [{"molecule_chembl_id": "CHEMBL1"}], "page_meta": {"next": "/chembl/api/data/molecule.json?offset=1"}})
            return Response({"molecules": [{"molecule_chembl_id": "CHEMBL2"}], "page_meta": {"next": None}})

    session = Session()
    records = fetch_approved_molecules(session)
    assert [record["molecule_chembl_id"] for record in records] == ["CHEMBL1", "CHEMBL2"]
    assert session.urls[1][0] == "https://www.ebi.ac.uk/chembl/api/data/molecule.json?offset=1"
    assert session.urls[1][1] is None


def test_record_from_chembl_extracts_source_fields_and_cross_references():
    raw = {
        "molecule_chembl_id": "CHEMBL1",
        "pref_name": "Example",
        "molecule_type": "Small molecule",
        "max_phase": "4.0",
        "molecule_structures": {"canonical_smiles": "CCO"},
        "molecule_hierarchy": {"parent_chembl_id": "CHEMBL1"},
        "cross_references": [
            {"xref_src": "DrugBank", "xref_id": "DB00001"},
            {"xref_src": "PubChem", "xref_id": "CID: 702"},
        ],
    }
    record = record_from_chembl(raw)
    assert record is not None
    assert record.chembl_id == "CHEMBL1"
    assert record.max_phase == 4
    assert record.drugbank_id == "DB00001"
    assert record.pubchem_cid == "702"


def test_normalise_and_deduplicate_rejects_bad_smiles_large_molecules_and_salts():
    records = [
        LigandRecord("CHEMBL1", "Parent", "CCO", "Small molecule", 4, "CHEMBL1"),
        LigandRecord("CHEMBL2", "Salt", "CCO.Cl", "Small molecule", 4, "CHEMBL1"),
        LigandRecord("CHEMBL3", "Broken", "not smiles", "Small molecule", 4, "CHEMBL3"),
        LigandRecord("CHEMBL4", "Huge", "C" * 200, "Small molecule", 4, "CHEMBL4"),
    ]
    results = {record.chembl_id: record for record in normalise_and_deduplicate(records, max_molecular_weight=100)}
    assert results["CHEMBL1"].prep_status == "pending"
    assert results["CHEMBL1"].prepared_smiles == "CCO"
    assert results["CHEMBL2"].prep_status == "excluded_duplicate"
    assert results["CHEMBL3"].prep_status == "rejected_unparseable"
    assert results["CHEMBL4"].prep_status == "rejected_size"


def test_normalise_resolves_known_chembl_parent_and_flags_unresolved_multifragment_records():
    records = [
        LigandRecord("CHEMBL_PARENT", "Parent", "CCN", "Small molecule", 4, "CHEMBL_PARENT"),
        LigandRecord("CHEMBL_SALT", "Parent salt", "CCN.O=S(=O)(O)c1ccccc1", "Small molecule", 4, "CHEMBL_PARENT"),
        LigandRecord(
            "CHEMBL_MIX",
            "Unresolved mix",
            "CCN.O=C(O)c1ccc2ccccc2c1O",
            "Small molecule",
            4,
            "CHEMBL_MIX",
        ),
    ]
    results = {record.chembl_id: record for record in normalise_and_deduplicate(records)}
    assert results["CHEMBL_PARENT"].prep_status == "pending"
    assert results["CHEMBL_SALT"].prep_status == "excluded_duplicate"
    assert "CHEMBL_PARENT" in results["CHEMBL_SALT"].prep_notes
    assert results["CHEMBL_MIX"].prep_status == "rejected_unresolved_multifragment"


def test_normalise_rejects_carbon_free_structures_as_inorganic():
    records = [
        LigandRecord("CHEMBL_INORGANIC", "Sodium chloride", "[Na+].[Cl-]", "Small molecule", 4, "CHEMBL_INORGANIC")
    ]
    result = normalise_and_deduplicate(records)[0]
    assert result.prep_status == "rejected_inorganic"


def test_prepare_library_reuses_preparation_and_manifest_tracks_success(monkeypatch, tmp_path):
    def fake_sdf(smiles, name, out_path, seed):
        assert smiles == "CCO"
        assert name == "CHEMBL1"
        assert seed == 42
        out_path.write_text("fake sdf")

    def fake_pdbqt(sdf_path, out_path):
        assert sdf_path.exists()
        out_path.write_text("ROOT\nATOM\nENDROOT\n")

    monkeypatch.setattr("molecular_docking.build_ligand_library.build_ligand_sdf", fake_sdf)
    monkeypatch.setattr("molecular_docking.build_ligand_library.prepare_ligand_pdbqt", fake_pdbqt)
    record = LigandRecord(
        "CHEMBL1", "Example", "CCO", "Small molecule", 4, "CHEMBL1", prepared_smiles="CCO", molecular_weight=46.07
    )
    completed = prepare_library([record], tmp_path)
    assert completed[0].prep_status == "ok"
    assert (tmp_path / completed[0].pdbqt_path).exists()
    assert not (tmp_path / ".staging_sdf").exists()
    manifest = write_manifest(completed, tmp_path)
    assert "CHEMBL1" in manifest.read_text()
    assert Path(completed[0].pdbqt_path).suffix == ".pdbqt"
