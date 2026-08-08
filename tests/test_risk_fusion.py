import pytest

from risk_fusion.fuse_scores import (
    check_mdr,
    docking_component,
    fuse,
    genomic_component,
    score_band,
)


def _genomic_result(**overrides) -> dict:
    defaults = dict(
        gene="katG",
        drug="Isoniazid",
        predicted_phenotype="resistant",
        driving_mutation="S315T",
        who_confidence_grade=1,
        confidence_label="Assoc w R",
        resistant_mutations_found=["S315T"],
        positions_checked=5,
        calls=[
            {
                "mutation": "S315T",
                "observed_aa": "T",
                "status": "resistant",
                "who_confidence_grade": 1,
                "confidence_label": "Assoc w R",
            },
            {
                "mutation": "W328L",
                "observed_aa": "A",
                "status": "novel_variant",
                "who_confidence_grade": 3,
                "confidence_label": "Uncertain significance",
            },
        ],
    )
    defaults.update(overrides)
    return defaults


def _susceptible_genomic_result(**overrides) -> dict:
    return _genomic_result(
        predicted_phenotype="susceptible",
        driving_mutation=None,
        who_confidence_grade=None,
        confidence_label=None,
        resistant_mutations_found=[],
        calls=[
            {
                "mutation": "S315T",
                "observed_aa": "S",
                "status": "wild_type",
                "who_confidence_grade": 1,
                "confidence_label": "Assoc w R",
            },
        ],
        **overrides,
    )


def _docking_result(**overrides) -> dict:
    defaults = dict(
        gene="katG",
        drug="Isoniazid",
        pubchem_cid=3767,
        pdb_id="2CCA",
        box_center=(84.4, 38.4, 49.1),
        box_size=(20.0, 20.0, 20.0),
        poses=[],
        best_affinity_kcal_per_mol=-5.8,
        docked_pdbqt_path="unused.pdbqt",
    )
    defaults.update(overrides)
    return defaults


def test_score_band_thresholds():
    assert score_band(95) == "High confidence resistant"
    assert score_band(85) == "High confidence resistant"
    assert score_band(84.9) == "Likely resistant"
    assert score_band(60) == "Likely resistant"
    assert score_band(35) == "Uncertain / insufficient evidence"
    assert score_band(0) == "Likely susceptible"


def test_genomic_component_grade1_resistant_scores_high():
    score, notes = genomic_component(_genomic_result())
    assert score == 95
    assert any("S315T" in n for n in notes)


def test_genomic_component_grade2_resistant_scores_lower_than_grade1():
    score, _ = genomic_component(_genomic_result(who_confidence_grade=2, confidence_label="Assoc w R - Interim"))
    assert score == 75


def test_genomic_component_susceptible_is_low():
    score, notes = genomic_component(_susceptible_genomic_result())
    assert score == 8
    assert notes == []


def test_genomic_component_susceptible_with_novel_variant_bumps_score():
    result = _genomic_result(
        predicted_phenotype="susceptible",
        driving_mutation=None,
        who_confidence_grade=None,
        confidence_label=None,
        resistant_mutations_found=[],
        calls=[
            {
                "mutation": "S315T",
                "observed_aa": "S",
                "status": "wild_type",
                "who_confidence_grade": 1,
                "confidence_label": "Assoc w R",
            },
            {
                "mutation": "W328L",
                "observed_aa": "Q",
                "status": "novel_variant",
                "who_confidence_grade": 3,
                "confidence_label": "Uncertain significance",
            },
        ],
    )
    score, notes = genomic_component(result)
    assert score == 8 + 7
    assert any("W328L" in n for n in notes)


def test_docking_component_returns_zero_when_no_docking_result():
    delta, notes = docking_component(None)
    assert delta == 0.0
    assert notes == []


def test_docking_component_weaker_than_reference_increases_score():
    delta, notes = docking_component(_docking_result(best_affinity_kcal_per_mol=-4.0))
    assert delta > 0
    assert notes


def test_docking_component_stronger_than_reference_decreases_score():
    delta, _ = docking_component(_docking_result(best_affinity_kcal_per_mol=-10.0))
    assert delta < 0


def test_docking_component_is_capped():
    delta, _ = docking_component(_docking_result(best_affinity_kcal_per_mol=0.0))
    assert delta == pytest.approx(8.0)
    delta2, _ = docking_component(_docking_result(best_affinity_kcal_per_mol=-30.0))
    assert delta2 == pytest.approx(-8.0)


def test_fuse_resistant_result_has_high_score_and_alternatives():
    result = fuse(_genomic_result())
    assert result.band in ("High confidence resistant", "Likely resistant")
    assert "Rifampicin" in result.alternative_drugs


def test_fuse_susceptible_result_has_no_alternatives():
    result = fuse(_susceptible_genomic_result())
    assert result.band == "Likely susceptible"
    assert result.alternative_drugs == []


def test_fuse_mutation_report_flags_driving_mutation():
    result = fuse(_genomic_result())
    driving_rows = [r for r in result.mutation_report if r["is_driving_mutation"]]
    assert len(driving_rows) == 1
    assert driving_rows[0]["mutation"] == "S315T"


def test_fuse_docking_cannot_flip_susceptible_to_resistant():
    # Docking's max adjustment (+8) must never be enough by itself to cross
    # into a resistant band from the susceptible base score.
    result = fuse(_susceptible_genomic_result(), _docking_result(best_affinity_kcal_per_mol=100.0))
    assert result.band in ("Likely susceptible", "Uncertain / insufficient evidence")


def test_fuse_attaches_public_health_context_without_affecting_score():
    ctx = {
        "country": "Zambia",
        "year": 2023,
        "population": 20_723_959,
        "incident_cases": 59000,
        "case_fatality_ratio": 0.09,
        "treatment_success_rate_pct": 92.0,
    }
    with_ctx = fuse(_genomic_result(), public_health_context=ctx)
    without_ctx = fuse(_genomic_result())
    assert with_ctx.resistance_score == without_ctx.resistance_score
    assert with_ctx.public_health_context == ctx


def test_check_mdr_flags_when_both_resistant():
    katg = fuse(_genomic_result(gene="katG", drug="Isoniazid"))
    rpob = fuse(_genomic_result(gene="rpoB", drug="Rifampicin", driving_mutation="S450L"))
    mdr = check_mdr(katg, rpob)
    assert mdr is not None
    assert mdr["mdr_flag"] is True
    assert "Bedaquiline" in mdr["recommended_regimen"]


def test_check_mdr_returns_none_when_only_one_resistant():
    katg = fuse(_genomic_result(gene="katG", drug="Isoniazid"))
    rpob = fuse(_susceptible_genomic_result(gene="rpoB", drug="Rifampicin"))
    assert check_mdr(katg, rpob) is None
