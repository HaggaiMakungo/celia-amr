"""Unified command-line workflow for the Celia AMR prototype.

The individual modules remain runnable on their own.  This module is the
single entry point that connects them for a sample-level assessment:

    celia analyze --gene katG --fasta genomic_analysis/examples/katG_resistant_example.fasta

Docking is deliberately opt-in because it can download structures and run a
native binary.  Zambia's WHO burden data is also opt-in: it is useful context,
but never changes the per-sample resistance score.
"""

import argparse
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Callable

from genomic_analysis.analyze_mutations import analyze as analyze_genome
from risk_fusion.fuse_scores import FusionResult, fuse


def public_health_context_from_snapshot(snapshot: Any) -> dict:
    """Convert a spread-model burden snapshot into reportable context."""
    return {
        "country": snapshot.country,
        "year": snapshot.year,
        "population": snapshot.population,
        "incident_cases": snapshot.incident_cases,
        "case_fatality_ratio": snapshot.case_fatality_ratio,
        "treatment_success_rate_pct": snapshot.treatment_success_rate_pct,
    }


def run_analysis(
    gene: str,
    fasta_path: Path,
    *,
    with_docking: bool = False,
    with_burden_context: bool = False,
    burden_year: int = 2023,
    docking_runner: Callable[[str], dict] | None = None,
) -> tuple[FusionResult, dict, dict | None]:
    """Run the connected AMR workflow and return its three evidence layers.

    The return value is ``(assessment, genomic_result, docking_result)``.
    ``assessment.public_health_context`` contains the optional WHO context.
    ``docking_runner`` is injectable for programmatic callers and tests; normal
    CLI runs load the real AutoDock Vina workflow only when requested.
    """
    genomic_result = analyze_genome(gene, fasta_path)

    docking_result = None
    if with_docking:
        if docking_runner is None:
            # Keep these heavy/native dependencies out of normal genome-only
            # runs and make failures local to the requested optional feature.
            from molecular_docking.run_docking import dock

            docking_runner = dock
        docking_result = docking_runner(gene)

    public_health_context = None
    if with_burden_context:
        from spread_modelling.seir_model import load_burden_snapshot

        public_health_context = public_health_context_from_snapshot(
            load_burden_snapshot(burden_year)
        )

    assessment = fuse(genomic_result, docking_result, public_health_context)
    return assessment, genomic_result, docking_result


def report_as_dict(
    assessment: FusionResult,
    genomic_result: dict,
    docking_result: dict | None,
) -> dict:
    """Build a JSON-safe report while keeping every evidence source visible."""
    return {
        "assessment": asdict(assessment),
        "genomic_analysis": genomic_result,
        "molecular_docking": _json_safe(docking_result),
        "public_health_context": assessment.public_health_context,
    }


def _json_safe(value: Any) -> Any:
    if is_dataclass(value):
        return _json_safe(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def print_text_report(assessment: FusionResult, *, docking_requested: bool) -> None:
    """Print a concise, readable assessment for terminal users."""
    print("Celia AMR assessment")
    print(f"Gene: {assessment.gene}  Drug: {assessment.drug}")
    print(f"Resistance score: {assessment.resistance_score:.1f}/100 ({assessment.band})")
    print(
        "Evidence: "
        f"genomic {assessment.genomic_component:.1f} points; "
        f"docking {assessment.docking_component:+.1f} points"
    )
    for note in assessment.notes:
        print(f"  - {note}")

    print("Mutation report:")
    for row in assessment.mutation_report:
        driving = " [driving mutation]" if row["is_driving_mutation"] else ""
        print(
            f"  - {row['mutation']}: {row['status']} "
            f"(WHO grade {row['who_confidence_grade']}){driving}"
        )

    if assessment.alternative_drugs:
        print("Alternative drugs to consider: " + ", ".join(assessment.alternative_drugs))
    if not docking_requested:
        print("Docking: not run (add --with-docking to include it).")
    if assessment.public_health_context:
        context = assessment.public_health_context
        print(
            f"Public-health context: {context['country']} ({context['year']}), "
            f"{context['incident_cases']:,.0f} estimated incident TB cases, "
            f"CFR {context['case_fatality_ratio']:.1%}."
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="celia",
        description="Celia AMR: integrated tuberculosis resistance assessment.",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    analyze = subcommands.add_parser(
        "analyze",
        help="combine genomic evidence with optional docking and WHO context",
    )
    analyze.add_argument("--gene", required=True, choices=["katG", "rpoB"])
    analyze.add_argument("--fasta", required=True, type=Path, help="protein FASTA for the selected gene")
    analyze.add_argument(
        "--with-docking",
        action="store_true",
        help="run real AutoDock Vina docking (may download data and take longer)",
    )
    analyze.add_argument(
        "--with-burden-context",
        action="store_true",
        help="attach Zambia WHO TB burden context; it does not change the resistance score",
    )
    analyze.add_argument("--burden-year", type=int, default=2023)
    analyze.add_argument("--format", choices=["text", "json"], default="text")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command != "analyze":  # Defensive guard if future subcommands are added.
        raise ValueError(f"Unsupported command: {args.command}")

    assessment, genomic_result, docking_result = run_analysis(
        args.gene,
        args.fasta,
        with_docking=args.with_docking,
        with_burden_context=args.with_burden_context,
        burden_year=args.burden_year,
    )
    if args.format == "json":
        print(json.dumps(report_as_dict(assessment, genomic_result, docking_result), indent=2))
    else:
        print_text_report(assessment, docking_requested=args.with_docking)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
