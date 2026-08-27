from molecular_docking.rank_ligands import ranked_rows, review_flags


def _row(chembl_id, affinity, smiles, status="ok"):
    return {
        "chembl_id": chembl_id, "pref_name": chembl_id, "prep_status": "ok",
        "docking_status": status, "best_affinity_kcal_per_mol": str(affinity), "prepared_smiles": smiles,
    }


def test_ranked_rows_orders_more_negative_affinities_first_and_skips_failures():
    ranked = ranked_rows([_row("CHEMBL2", -7.0, "CCO"), _row("CHEMBL1", -9.0, "CCN"), _row("CHEMBL3", -12.0, "CC", "failed")], 1)
    assert [row["chembl_id"] for row in ranked] == ["CHEMBL1", "CHEMBL2"]
    assert ranked[0]["second_pass_recommendation"] == "rerun_high_exhaustiveness"
    assert ranked[1]["second_pass_recommendation"] == "not_shortlisted"


def test_review_flags_are_context_not_exclusions():
    flags = review_flags({"molecular_weight": 700.0, "rdkit_logp": 6.0, "tpsa": 150.0, "hbd": 1, "hba": 2, "rotatable_bonds": 13, "formal_charge": 2})
    assert flags == ["high_molecular_weight", "high_flexibility", "high_formal_charge", "high_logp", "high_polar_surface_area"]
