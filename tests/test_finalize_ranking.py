from molecular_docking.finalize_ranking import final_rows


def test_final_rows_replaces_only_confirmed_second_pass_scores_and_reranks():
    first = [
        {"rank": "1", "chembl_id": "CHEMBL1", "pref_name": "One", "best_affinity_kcal_per_mol": "-10.0"},
        {"rank": "2", "chembl_id": "CHEMBL2", "pref_name": "Two", "best_affinity_kcal_per_mol": "-9.0"},
    ]
    second = [
        {"chembl_id": "CHEMBL2", "second_pass_status": "ok", "second_pass_best_affinity_kcal_per_mol": "-11.0"},
        {"chembl_id": "CHEMBL1", "second_pass_status": "failed", "second_pass_best_affinity_kcal_per_mol": "-12.0"},
    ]
    ranked = final_rows(first, second)
    assert [row["chembl_id"] for row in ranked] == ["CHEMBL2", "CHEMBL1"]
    assert ranked[0]["affinity_source"] == "second_pass_exhaustiveness_8"
    assert ranked[1]["affinity_source"] == "first_pass_exhaustiveness_4"
