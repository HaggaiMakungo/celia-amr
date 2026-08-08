import pytest

from spread_modelling.seir_model import (
    BURDEN_CSV_PATH,
    OUTCOMES_CSV_PATH,
    BurdenSnapshot,
    build_parameters,
    calibrate_beta,
    load_burden_snapshot,
    seir_rhs,
    simulate,
)


def test_burden_data_files_exist():
    assert BURDEN_CSV_PATH.exists()
    assert OUTCOMES_CSV_PATH.exists()


def test_load_burden_snapshot_reads_real_2023_data():
    snapshot = load_burden_snapshot(2023)
    assert snapshot.country == "Zambia"
    assert snapshot.year == 2023
    # Sanity-check against the actual WHO-reported magnitude, not exact
    # values, so this test doesn't silently rot if WHO revises estimates.
    assert 15_000_000 < snapshot.population < 30_000_000
    assert 30_000 < snapshot.incident_cases < 100_000
    assert 0 < snapshot.case_fatality_ratio < 1


def test_load_burden_snapshot_rejects_unavailable_year():
    with pytest.raises(ValueError):
        load_burden_snapshot(1800)


def test_calibrate_beta_solves_flow_balance():
    beta = calibrate_beta(N=1000.0, S0=800.0, E0=100.0, I0=50.0, sigma=0.01)
    # beta * S0 * I0 / N should reproduce sigma * E0
    assert beta * 800.0 * 50.0 / 1000.0 == pytest.approx(0.01 * 100.0)


def test_calibrate_beta_rejects_non_positive_compartments():
    with pytest.raises(ValueError):
        calibrate_beta(N=1000.0, S0=0.0, E0=100.0, I0=50.0, sigma=0.01)


def _fake_snapshot(**overrides) -> BurdenSnapshot:
    defaults = dict(
        country="Testland",
        year=2023,
        population=1_000_000.0,
        incident_cases=1_000.0,
        case_fatality_ratio=0.1,
        treatment_success_rate_pct=90.0,
    )
    defaults.update(overrides)
    return BurdenSnapshot(**defaults)


def test_build_parameters_matches_reported_incidence_at_t0():
    snapshot = _fake_snapshot()
    params = build_parameters(snapshot, sigma=0.01, disease_duration_years=1.0)
    # sigma * E0 (inflow to I) should equal the real reported incident cases.
    assert params.sigma * params.E0 == pytest.approx(snapshot.incident_cases)
    # (gamma + mu_tb) * I0 (outflow from I) should also equal it.
    assert (params.gamma + params.mu_tb) * params.I0 == pytest.approx(snapshot.incident_cases)


def test_build_parameters_splits_exit_rate_by_case_fatality_ratio():
    snapshot = _fake_snapshot(case_fatality_ratio=0.2)
    params = build_parameters(snapshot, sigma=0.01, disease_duration_years=1.0)
    exit_rate = params.gamma + params.mu_tb
    assert params.mu_tb == pytest.approx(0.2 * exit_rate)
    assert params.gamma == pytest.approx(0.8 * exit_rate)


def test_build_parameters_rejects_oversized_latent_pool():
    # An implausibly tiny sigma inflates E0 past the population size.
    snapshot = _fake_snapshot(population=1_000.0, incident_cases=1_000.0)
    with pytest.raises(ValueError):
        build_parameters(snapshot, sigma=0.0001, disease_duration_years=1.0)


def test_seir_rhs_population_declines_only_by_tb_death_flow():
    # No births/other-cause deaths are modeled, so the only way S+E+I+R
    # changes is TB deaths (mu_tb * I) leaving the living population.
    snapshot = _fake_snapshot()
    params = build_parameters(snapshot)
    y = [params.S0, params.E0, params.I0, params.R0]
    derivatives = seir_rhs(0.0, y, params)
    assert sum(derivatives) == pytest.approx(-params.mu_tb * params.I0)


def test_simulate_produces_monotonic_time_and_starts_at_initial_conditions():
    snapshot = _fake_snapshot()
    params = build_parameters(snapshot)
    trajectory = simulate(params, years=5.0, points_per_year=4)

    assert trajectory["t_years"].is_monotonic_increasing
    assert trajectory.iloc[0]["S"] == pytest.approx(params.S0)
    assert trajectory.iloc[0]["I"] == pytest.approx(params.I0)


def test_simulate_population_declines_due_to_tb_deaths():
    snapshot = _fake_snapshot()
    params = build_parameters(snapshot)
    trajectory = simulate(params, years=5.0, points_per_year=4)

    totals = trajectory["S"] + trajectory["E"] + trajectory["I"] + trajectory["R"]
    # Monotonically non-increasing: the only outflow from the living
    # population is TB death, there's no inflow to offset it.
    assert (totals.diff().dropna() <= 1e-9).all()
    assert totals.iloc[-1] < totals.iloc[0]
