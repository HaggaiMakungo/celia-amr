import pytest
import subprocess

from molecular_docking.run_docking import (
    TARGETS,
    DockingPose,
    clean_receptor_pdb,
    compute_box_center,
    parse_vina_output,
    run_vina,
)

VINA_STDOUT_SAMPLE = """\
Scoring function : vina
Rigid receptor: 2CCA_receptor.pdbqt
Ligand: isoniazid.pdbqt
Grid center: X 84.403 Y 38.356 Z 49.07
Grid size  : X 20 Y 20 Z 20
Grid space : 0.375
Exhaustiveness: 8
CPU: 0
Verbosity: 1

Computing Vina grid ... done.
Performing docking (random seed: 42) ...
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -5.858          0          0
   2       -5.782       7.22      9.063
   3       -5.714      1.705      2.025
"""


def test_targets_registry_covers_both_genes():
    assert set(TARGETS.keys()) == {"katG", "rpoB"}
    for gene, target in TARGETS.items():
        assert target.gene == gene
        assert target.pdb_id
        assert target.drug_name


def test_parse_vina_output_extracts_all_poses():
    poses = parse_vina_output(VINA_STDOUT_SAMPLE)
    assert len(poses) == 3
    assert poses[0] == DockingPose(mode=1, affinity_kcal_per_mol=-5.858, rmsd_lower_bound=0.0, rmsd_upper_bound=0.0)
    assert poses[1].mode == 2
    assert poses[1].affinity_kcal_per_mol == pytest.approx(-5.782)
    assert poses[1].rmsd_lower_bound == pytest.approx(7.22)
    assert poses[1].rmsd_upper_bound == pytest.approx(9.063)


def test_parse_vina_output_ignores_non_pose_lines():
    poses = parse_vina_output("some\nrandom\ntext\nwith no pose rows")
    assert poses == []


def test_run_vina_reports_a_clear_timeout(monkeypatch, tmp_path):
    def timeout(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd="vina", timeout=15)

    monkeypatch.setattr("molecular_docking.run_docking.subprocess.run", timeout)
    with pytest.raises(RuntimeError, match="timed out after 15 seconds"):
        run_vina(
            tmp_path / "vina.exe",
            tmp_path / "receptor.pdbqt",
            tmp_path / "ligand.pdbqt",
            (0.0, 0.0, 0.0),
            (20.0, 20.0, 20.0),
            tmp_path / "result.pdbqt",
            timeout_seconds=15,
        )


def test_parse_vina_output_poses_are_sorted_best_first_by_construction():
    poses = parse_vina_output(VINA_STDOUT_SAMPLE)
    affinities = [p.affinity_kcal_per_mol for p in poses]
    assert affinities == sorted(affinities)  # most negative (best) first


def test_clean_receptor_pdb_strips_only_requested_hetatm_residues(tmp_path):
    pdb_path = tmp_path / "sample.pdb"
    pdb_path.write_text(
        "ATOM      1  N   ALA A   1      11.000  12.000  13.000  1.00  0.00           N\n"
        "HETATM 1000  O   HOH A2001      59.372   5.076  50.482  1.00 39.62           O\n"
        "HETATM 1001  FE  HEM A1741      84.406  37.968  51.393  1.00 16.46          FE\n"
        "END\n"
    )
    out_path = tmp_path / "cleaned.pdb"
    clean_receptor_pdb(pdb_path, strip_resnames=("HOH",), out_path=out_path)

    cleaned = out_path.read_text()
    assert "ALA" in cleaned
    assert "HEM" in cleaned
    assert "HOH" not in cleaned


def test_compute_box_center_averages_real_coordinates(tmp_path):
    pdb_path = tmp_path / "sample.pdb"
    pdb_path.write_text(
        "HETATM 1000  FE  HEM A1741      10.000  20.000  30.000  1.00 16.46          FE\n"
        "HETATM 1001  NA  HEM A1741      20.000  30.000  40.000  1.00 16.46           N\n"
        "HETATM 1002  NA  HEM B1741       0.000   0.000   0.000  1.00 16.46           N\n"
    )
    center = compute_box_center(pdb_path, resname="HEM", chain="A")
    assert center == pytest.approx((15.0, 25.0, 35.0))


def test_compute_box_center_raises_when_no_matching_atoms(tmp_path):
    pdb_path = tmp_path / "sample.pdb"
    pdb_path.write_text("ATOM      1  N   ALA A   1      11.000  12.000  13.000  1.00  0.00           N\n")
    with pytest.raises(ValueError):
        compute_box_center(pdb_path, resname="HEM", chain="A")
