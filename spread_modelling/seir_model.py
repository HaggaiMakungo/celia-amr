"""Branch 2: Spread Modelling.

A compartmental SEIR (Susceptible - Exposed/latent - Infectious - Recovered)
model for TB transmission, parameterized from real WHO TB burden data for
Zambia (`data/tb_burden/zambia_tb_burden.csv`, `zambia_tb_outcomes.csv`).

Standalone script — no dependency on the other three branches. Run directly:

    python seir_model.py --year 2023 --project-years 10

## Model

    dS/dt = -beta * S * I / N
    dE/dt =  beta * S * I / N - sigma * E
    dI/dt =  sigma * E - (gamma + mu_tb) * I
    dR/dt =  gamma * I

S: susceptible, E: latently infected (not infectious), I: active infectious
disease, R: recovered/treated (assumed immune to reinfection in this
first-pass model). `N` is fixed at the reference year's population and used
only as the transmission-mixing denominator (a standard frequency-dependent
transmission simplification). It is *not* a conserved total: `S+E+I+R`
declines over time as `mu_tb * I` removes TB deaths from the living
population, and this prototype does not model births, non-TB deaths, or
migration replacing them.

## Where each parameter comes from

- `N` (population) and `e_inc_num` (estimated incident active TB cases in the
  reference year) are read directly from the WHO burden CSV for Zambia.
- `cfr` (case fatality ratio: the fraction of incident cases that die) is
  read directly from the same CSV. It splits the exit rate from the
  infectious compartment `I` into a death rate `mu_tb` and a recovery rate
  `gamma`.
- `disease_duration_years` (average time an individual spends in the
  infectious compartment before death, cure, or treatment completion) is a
  literature-standard assumption (~1 year, combining diagnostic delay and
  treatment course; see Vynnycky & Fine 1997, Blower et al. 1995 for typical
  ranges of 0.5-2 years in similar models) — WHO burden data does not report
  disease duration directly, so this is the one parameter here that is a
  modeling assumption rather than a direct data read.
- `sigma` (annual per-latent-individual progression rate from `E` to `I`) is
  likewise a literature-informed assumption, not a WHO-reported figure.
  General-population lifetime progression risk is often cited near 5-10%,
  but that figure is dominated by decades of slow reactivation risk in
  low-transmission settings. In a high-transmission, high-recent-infection
  setting like Zambia, active disease is disproportionately driven by recent
  (fast) progression, which is better approximated by a higher annual rate;
  this model uses `sigma = 0.01/year` as a single-pathway approximation
  (see Blower et al. 1995 for fast/slow progression modeling). This is the
  most uncertain parameter in the model and is exposed as a CLI flag so it
  can be swept/sensitivity-tested.
- `beta` (the transmission coefficient) is not looked up anywhere — it is
  *calibrated* by this script so that the model's endemic incidence flow
  (`sigma * E0`) matches the real WHO-reported incident case count in the
  reference year, holding the other parameters fixed. See `calibrate_beta`.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

BURDEN_CSV_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "tb_burden" / "zambia_tb_burden.csv"
)
OUTCOMES_CSV_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "tb_burden" / "zambia_tb_outcomes.csv"
)

DEFAULT_SIGMA_PER_YEAR = 0.01
DEFAULT_DISEASE_DURATION_YEARS = 1.0


@dataclass(frozen=True)
class BurdenSnapshot:
    country: str
    year: int
    population: float
    incident_cases: float
    case_fatality_ratio: float
    treatment_success_rate_pct: float | None


def load_burden_snapshot(
    year: int,
    burden_csv: Path = BURDEN_CSV_PATH,
    outcomes_csv: Path = OUTCOMES_CSV_PATH,
) -> BurdenSnapshot:
    """Read the real WHO-reported burden figures for Zambia in `year`."""
    burden = pd.read_csv(burden_csv)
    row = burden[burden["year"] == year]
    if row.empty:
        available = sorted(burden["year"].tolist())
        raise ValueError(f"No burden data for year {year}. Available years: {available}")
    row = row.iloc[0]

    tsr = None
    try:
        outcomes = pd.read_csv(outcomes_csv)
        # Treatment outcomes for a given notification year are reported in
        # the *following* year's cohort file (cohorts are only closed out
        # after the full treatment course completes).
        outcome_row = outcomes[outcomes["year"] == year + 1]
        if not outcome_row.empty and pd.notna(outcome_row.iloc[0]["c_new_tsr"]):
            tsr = float(outcome_row.iloc[0]["c_new_tsr"])
    except FileNotFoundError:
        pass

    return BurdenSnapshot(
        country=row["country"],
        year=int(row["year"]),
        population=float(row["e_pop_num"]),
        incident_cases=float(row["e_inc_num"]),
        case_fatality_ratio=float(row["cfr"]),
        treatment_success_rate_pct=tsr,
    )


@dataclass(frozen=True)
class SEIRParameters:
    N: float
    beta: float
    sigma: float
    gamma: float
    mu_tb: float
    S0: float
    E0: float
    I0: float
    R0: float


def calibrate_beta(N: float, S0: float, E0: float, I0: float, sigma: float) -> float:
    """Solve beta * S0 * I0 / N = sigma * E0 for beta.

    This calibrates transmission so the model's initial new-infection flow
    into E matches the initial progression flow out of E (i.e. E starts at
    quasi-steady-state given the other fixed parameters), which in turn was
    sized to match the real reported incidence flow (see `build_parameters`).
    """
    if S0 <= 0 or I0 <= 0:
        raise ValueError("S0 and I0 must be positive to calibrate beta")
    return sigma * E0 * N / (S0 * I0)


def build_parameters(
    snapshot: BurdenSnapshot,
    sigma: float = DEFAULT_SIGMA_PER_YEAR,
    disease_duration_years: float = DEFAULT_DISEASE_DURATION_YEARS,
) -> SEIRParameters:
    """Derive a self-consistent SEIRParameters set from a real burden snapshot.

    E0 and I0 are sized so that, at t=0, the model's internal flows
    (sigma * E0 into I, and (gamma + mu_tb) * I0 out of I) both equal the
    real reported annual incident case count -- i.e. the model starts at the
    flow-balance point implied by the actual data, rather than at an
    arbitrary guess.
    """
    N = snapshot.population
    incident_cases = snapshot.incident_cases
    cfr = snapshot.case_fatality_ratio

    exit_rate = 1.0 / disease_duration_years  # gamma + mu_tb
    mu_tb = cfr * exit_rate
    gamma = (1 - cfr) * exit_rate

    E0 = incident_cases / sigma
    I0 = incident_cases / exit_rate
    R0 = 0.0
    S0 = N - E0 - I0 - R0
    if S0 <= 0:
        raise ValueError(
            f"Derived S0 <= 0 (E0={E0:.0f}, I0={I0:.0f} exceed population {N:.0f}); "
            "try a larger --sigma"
        )

    beta = calibrate_beta(N, S0, E0, I0, sigma)

    return SEIRParameters(N=N, beta=beta, sigma=sigma, gamma=gamma, mu_tb=mu_tb, S0=S0, E0=E0, I0=I0, R0=R0)


def seir_rhs(t: float, y: np.ndarray, params: SEIRParameters) -> list:
    S, E, I, R = y
    N = params.N
    dS = -params.beta * S * I / N
    dE = params.beta * S * I / N - params.sigma * E
    dI = params.sigma * E - (params.gamma + params.mu_tb) * I
    dR = params.gamma * I
    return [dS, dE, dI, dR]


def simulate(params: SEIRParameters, years: float, points_per_year: int = 12) -> pd.DataFrame:
    """Integrate the SEIR system forward from t=0 (the reference year) for `years`."""
    t_span = (0.0, years)
    t_eval = np.linspace(0.0, years, int(years * points_per_year) + 1)
    y0 = [params.S0, params.E0, params.I0, params.R0]

    solution = solve_ivp(seir_rhs, t_span, y0, args=(params,), t_eval=t_eval, method="RK45")
    if not solution.success:
        raise RuntimeError(f"SEIR integration failed: {solution.message}")

    return pd.DataFrame(
        {
            "t_years": solution.t,
            "S": solution.y[0],
            "E": solution.y[1],
            "I": solution.y[2],
            "R": solution.y[3],
            "incidence_flow": params.sigma * solution.y[1],
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="TB SEIR spread model (Zambia, WHO burden data)")
    parser.add_argument("--year", type=int, default=2023, help="Reference year for calibration")
    parser.add_argument("--project-years", type=float, default=10.0, help="Years to simulate forward")
    parser.add_argument("--sigma", type=float, default=DEFAULT_SIGMA_PER_YEAR, help="E->I progression rate/year")
    parser.add_argument(
        "--disease-duration-years",
        type=float,
        default=DEFAULT_DISEASE_DURATION_YEARS,
        help="Average time in the infectious compartment before exit (death/recovery)",
    )
    parser.add_argument("--output", type=Path, default=None, help="Optional CSV path to write the trajectory")
    args = parser.parse_args()

    snapshot = load_burden_snapshot(args.year)
    params = build_parameters(snapshot, sigma=args.sigma, disease_duration_years=args.disease_duration_years)
    trajectory = simulate(params, years=args.project_years)

    print(f"{snapshot.country} TB SEIR model, calibrated to {snapshot.year} WHO burden data")
    print(f"  Population: {snapshot.population:,.0f}")
    print(f"  Reported incident cases ({snapshot.year}): {snapshot.incident_cases:,.0f}")
    print(f"  Case fatality ratio: {snapshot.case_fatality_ratio:.1%}")
    if snapshot.treatment_success_rate_pct is not None:
        print(f"  Treatment success rate ({snapshot.year + 1} cohort): {snapshot.treatment_success_rate_pct:.0f}%")
    print(f"  Calibrated beta: {params.beta:.4f}/year   sigma: {params.sigma:.4f}/year")
    print(f"  gamma: {params.gamma:.4f}/year   mu_tb: {params.mu_tb:.4f}/year")
    print(f"  Initial conditions: S0={params.S0:,.0f} E0={params.E0:,.0f} I0={params.I0:,.0f} R0={params.R0:,.0f}")
    print()
    print(trajectory.iloc[[0, -1]].to_string(index=False))

    if args.output:
        trajectory.to_csv(args.output, index=False)
        print(f"\nFull trajectory written to {args.output}")


if __name__ == "__main__":
    main()
