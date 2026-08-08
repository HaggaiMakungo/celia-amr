"""Branch 4: Risk Fusion Engine.

Combines the real output shapes produced by the other three branches into a
single resistance assessment: a resistance score, a mutation report, and
alternative-drug suggestions. No external dependency, no ML model — this is
plain, documented, weighted/rule-based Python, as scoped in the project brief.

Standalone script — imports the other branches only inside `main()`, and only
the ones actually needed for the requested run (spread-modelling context is
always pulled in; docking is opt-in via `--with-docking` since it needs the
native-toolchain environment set up for `molecular_docking`). The `fuse()`
function itself takes plain dicts shaped like `genomic_analysis.analyze()`'s
and `molecular_docking.dock()`'s return values, so it can be tested and used
without importing either module.

Run directly:

    python fuse_scores.py --gene katG --fasta ../genomic_analysis/examples/katG_resistant_example.fasta

## Design decisions

- **Spread-modelling data is deliberately excluded from the resistance
  score itself.** Folding a population-level burden/incidence figure into an
  individual sample's resistance probability without a proper Bayesian prior
  would conflate population and individual-level inference (an ecological
  fallacy) — the WHO burden data instead surfaces as a separate
  `public_health_context` section, informative for public-health-org users
  without corrupting the per-sample genomic/docking evidence.
- **Genomic evidence is the primary, clinically-grounded signal.** The WHO
  catalogue confidence grade drives the base score. Docking affinity is a
  secondary, explicitly low-confidence, capped adjustment — Vina scores from
  a single rigid-receptor run are not validated resistance predictors on
  their own, so they can nudge the score only a little, not flip a call.
"""

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

# --- genomic evidence -> base resistance score ---------------------------

GENOMIC_BASE_SCORE = {
    1: 95,  # WHO grade 1: Assoc w R
    2: 75,  # WHO grade 2: Assoc w R - Interim
}
SUSCEPTIBLE_BASE_SCORE = 8
NOVEL_VARIANT_UNCERTAINTY_BUMP = 7  # a variant at a known site that isn't in the catalogue warrants caution, not a resistance call

# --- docking evidence -> small capped adjustment --------------------------

# Informal midpoint of typical drug-like Vina affinities (-5 to -12 kcal/mol);
# not a validated clinical cutoff. See module docstring for why this
# contribution is deliberately small.
DOCKING_REFERENCE_AFFINITY_KCAL = -7.0
DOCKING_KCAL_TO_SCORE = 2.0
DOCKING_MAX_ADJUSTMENT = 8

# --- score bands -----------------------------------------------------------

SCORE_BANDS = [
    (85, "High confidence resistant"),
    (60, "Likely resistant"),
    (35, "Uncertain / insufficient evidence"),
    (0, "Likely susceptible"),
]

ALTERNATIVE_FIRST_LINE_DRUGS = {
    "Isoniazid": ["Rifampicin", "Ethambutol", "Pyrazinamide"],
    "Rifampicin": ["Isoniazid", "Ethambutol", "Pyrazinamide"],
}

# WHO's current standard shorter regimen for multidrug/rifampicin-resistant
# TB (BPaLM), per the WHO 2022 rapid communication and 2022 consolidated
# guidelines on drug-resistant TB treatment.
MDR_TB_REGIMEN = ["Bedaquiline", "Pretomanid", "Linezolid", "Moxifloxacin"]


def clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def score_band(score: float) -> str:
    for threshold, label in SCORE_BANDS:
        if score >= threshold:
            return label
    return SCORE_BANDS[-1][1]


@dataclass(frozen=True)
class FusionResult:
    gene: str
    drug: str
    resistance_score: float
    band: str
    genomic_component: float
    docking_component: float
    mutation_report: list
    alternative_drugs: list
    notes: list = field(default_factory=list)
    public_health_context: dict | None = None


def genomic_component(genomic_result: dict) -> tuple[float, list]:
    """Score contribution + notes from a genomic_analysis.analyze() result."""
    notes = []
    phenotype = genomic_result["predicted_phenotype"]
    calls = genomic_result.get("calls", [])

    if phenotype == "resistant":
        grade = genomic_result["who_confidence_grade"]
        score = GENOMIC_BASE_SCORE.get(grade, 60)
        notes.append(
            f"Driving mutation {genomic_result['driving_mutation']} "
            f"(WHO grade {grade}: {genomic_result['confidence_label']})"
        )
    else:
        score = SUSCEPTIBLE_BASE_SCORE
        novel_variants = [c["mutation"] for c in calls if c["status"] == "novel_variant"]
        if novel_variants:
            score += NOVEL_VARIANT_UNCERTAINTY_BUMP
            notes.append(
                f"Novel variant(s) at cataloged position(s) not matching any known "
                f"allele: {', '.join(novel_variants)} -- not a resistance call, but uncataloged"
            )

    return score, notes


def docking_component(docking_result: dict | None) -> tuple[float, list]:
    """Small, capped score adjustment from a molecular_docking.dock() result."""
    if docking_result is None:
        return 0.0, []

    affinity = docking_result["best_affinity_kcal_per_mol"]
    if affinity is None:
        return 0.0, ["Docking produced no poses; no docking adjustment applied"]

    # Weaker binding (less negative / numerically larger affinity than the
    # reference) should push the score up (more resistance-leaning); stronger
    # binding should push it down. That means the delta is
    # (affinity - reference), not (reference - affinity).
    raw_delta = (affinity - DOCKING_REFERENCE_AFFINITY_KCAL) * DOCKING_KCAL_TO_SCORE
    delta = clamp(raw_delta, -DOCKING_MAX_ADJUSTMENT, DOCKING_MAX_ADJUSTMENT)
    notes = [
        f"Docking best affinity {affinity:.2f} kcal/mol vs. reference "
        f"{DOCKING_REFERENCE_AFFINITY_KCAL:.1f} kcal/mol -> {delta:+.1f} point adjustment "
        "(low-confidence signal, not independently validated)"
    ]
    return delta, notes


def build_mutation_report(genomic_result: dict) -> list:
    report = []
    driving = genomic_result.get("driving_mutation")
    for call in genomic_result.get("calls", []):
        report.append(
            {
                "mutation": call["mutation"],
                "observed_aa": call["observed_aa"],
                "status": call["status"],
                "who_confidence_grade": call["who_confidence_grade"],
                "confidence_label": call["confidence_label"],
                "is_driving_mutation": call["mutation"] == driving,
            }
        )
    return report


def build_alternative_drugs(drug: str, band: str) -> list:
    if band in ("High confidence resistant", "Likely resistant"):
        return list(ALTERNATIVE_FIRST_LINE_DRUGS.get(drug, []))
    return []


def fuse(
    genomic_result: dict,
    docking_result: dict | None = None,
    public_health_context: dict | None = None,
) -> FusionResult:
    gene = genomic_result["gene"]
    drug = genomic_result["drug"]

    g_score, g_notes = genomic_component(genomic_result)
    d_score, d_notes = docking_component(docking_result)

    total = clamp(g_score + d_score, 0, 100)
    band = score_band(total)

    return FusionResult(
        gene=gene,
        drug=drug,
        resistance_score=total,
        band=band,
        genomic_component=g_score,
        docking_component=d_score,
        mutation_report=build_mutation_report(genomic_result),
        alternative_drugs=build_alternative_drugs(drug, band),
        notes=g_notes + d_notes,
        public_health_context=public_health_context,
    )


def check_mdr(katg_result: FusionResult, rpob_result: FusionResult) -> dict | None:
    """Flag MDR-TB (resistant to both first-line drugs) and the WHO regimen."""
    resistant_bands = ("High confidence resistant", "Likely resistant")
    if katg_result.band in resistant_bands and rpob_result.band in resistant_bands:
        return {
            "mdr_flag": True,
            "recommended_regimen": list(MDR_TB_REGIMEN),
            "note": "Resistant to both Isoniazid and Rifampicin (MDR-TB); "
            "first-line substitution is not appropriate -- WHO's standard "
            "shorter MDR/RR-TB regimen (BPaLM) is recommended instead.",
        }
    return None


def _public_health_context_from_snapshot(snapshot) -> dict:
    return {
        "country": snapshot.country,
        "year": snapshot.year,
        "population": snapshot.population,
        "incident_cases": snapshot.incident_cases,
        "case_fatality_ratio": snapshot.case_fatality_ratio,
        "treatment_success_rate_pct": snapshot.treatment_success_rate_pct,
    }


def _print_result(result: FusionResult) -> None:
    print(f"Gene: {result.gene}  Drug: {result.drug}")
    print(f"Resistance score: {result.resistance_score:.1f}/100  ({result.band})")
    print(f"  genomic component: {result.genomic_component:.1f}   docking component: {result.docking_component:+.1f}")
    for note in result.notes:
        print(f"  note: {note}")
    print("Mutation report:")
    for row in result.mutation_report:
        marker = " <- driving" if row["is_driving_mutation"] else ""
        print(f"  {row['mutation']}: {row['status']} (WHO grade {row['who_confidence_grade']}){marker}")
    if result.alternative_drugs:
        print(f"Alternative drugs to consider: {', '.join(result.alternative_drugs)}")
    if result.public_health_context:
        ctx = result.public_health_context
        print(
            f"Public health context ({ctx['country']}, {ctx['year']}): "
            f"{ctx['incident_cases']:,.0f} incident cases, "
            f"CFR {ctx['case_fatality_ratio']:.1%}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="TB risk fusion engine")
    parser.add_argument("--gene", required=True, choices=["katG", "rpoB"])
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--with-docking", action="store_true", help="Also run real molecular docking (slow, needs native toolchain)")
    parser.add_argument("--with-burden-context", action="store_true", help="Attach WHO Zambia TB burden context")
    parser.add_argument("--burden-year", type=int, default=2023)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from genomic_analysis.analyze_mutations import analyze as analyze_genomic

    genomic_result = analyze_genomic(args.gene, args.fasta)

    docking_result = None
    if args.with_docking:
        from molecular_docking.run_docking import dock

        docking_result = dock(args.gene)

    public_health_context = None
    if args.with_burden_context:
        from spread_modelling.seir_model import load_burden_snapshot

        snapshot = load_burden_snapshot(args.burden_year)
        public_health_context = _public_health_context_from_snapshot(snapshot)

    result = fuse(genomic_result, docking_result, public_health_context)
    _print_result(result)


if __name__ == "__main__":
    main()
