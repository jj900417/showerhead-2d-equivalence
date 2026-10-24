import csv
import math
from pathlib import Path

import pytest

from showerhead_equiv import Ring, annulus_bracket, equivalent_2d, solve_r_prime
from showerhead_equiv.cli import main, read_reactions, read_rings

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

# A synthetic ring used throughout the tests (not the geometry of any reactor in the article).
SYNTHETIC = dict(r=0.5, N=18, L=24.0, H=5.0)


def test_example_regression():
    """The synthetic example rings x the article's two reactions reproduce the reference outputs.

    expected_output.csv was generated with this package and cross-checked against an
    independent implementation of Eqs (3) and (4) (agreement to at least six significant figures).
    """
    rings = read_rings(EXAMPLES / "geometry_example.csv")
    reactions = read_reactions(EXAMPLES / "reactions_example.csv")
    with open(EXAMPLES / "expected_output.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    assert len(rings) == len(expected) == 6
    for ring, exp in zip(rings, expected):
        eq = equivalent_2d(ring)
        assert eq.r_prime == pytest.approx(float(exp["R_prime_mm"]), rel=1e-7)
        assert eq.area_ratio == pytest.approx(float(exp["wall_area_ratio_3D_to_2D"]), rel=1e-7)
        assert eq.gamma_eff(reactions[0][1]) == pytest.approx(float(exp["gamma_eff_R1"]), rel=1e-5)
        assert eq.gamma_eff(reactions[1][1]) == pytest.approx(float(exp["gamma_eff_R2"]), rel=1e-5)


def test_small_gap_asymptote():
    """Thin annulus limit: R' -> (3 N r^4 / 32 L)^(1/3), area ratio -> N r / 2L."""
    eq = equivalent_2d(Ring(r=0.5, N=3000, L=150, H=1))
    assert eq.r_prime == pytest.approx((3 * 3000 * 0.5**4 / (32 * 150)) ** (1 / 3), rel=1e-3)
    assert eq.r_prime == pytest.approx(0.48935834, rel=1e-7)
    assert eq.area_ratio == pytest.approx(3000 * 0.5 / (2 * 150))


def test_conductance_is_preserved():
    """At the returned R' the annulus carries exactly N r^4 (Eq. 3)."""
    eq = equivalent_2d(Ring(**SYNTHETIC))
    k = eq.ri / eq.ro
    annulus = eq.ro**4 - eq.ri**4 - (eq.ro**2 - eq.ri**2) ** 2 / math.log(1 / k)
    assert annulus == pytest.approx(SYNTHETIC["N"] * SYNTHETIC["r"] ** 4, rel=1e-9)


def test_scale_invariance():
    """Scaling r and L together scales R' alike and leaves k and the area ratio unchanged."""
    a = equivalent_2d(Ring(**SYNTHETIC))
    b = equivalent_2d(Ring(r=2 * SYNTHETIC["r"], N=SYNTHETIC["N"], L=2 * SYNTHETIC["L"], H=SYNTHETIC["H"]))
    assert b.r_prime == pytest.approx(2 * a.r_prime, rel=1e-9)
    assert b.k == pytest.approx(a.k, rel=1e-9)
    assert b.area_ratio == pytest.approx(a.area_ratio, rel=1e-12)


def test_series_fallback_is_continuous():
    k_switch = 1 - 5e-4
    exact = 1 - k_switch**4 - (1 - k_switch**2) ** 2 / math.log(1 / k_switch)
    assert annulus_bracket(k_switch + 1e-12) == pytest.approx(exact, rel=1e-6)


def test_residual_is_tiny():
    eq = equivalent_2d(Ring(**SYNTHETIC))
    assert abs(eq.residual_rel) < 1e-9


def test_no_equivalent_annulus():
    with pytest.raises(ValueError, match="no equivalent annulus"):
        solve_r_prime(N=10**9, r=1.0, L=2.0)


@pytest.mark.parametrize(
    "kwargs",
    [dict(r=0, N=4, L=10, H=1), dict(r=1, N=0, L=10, H=1), dict(r=11, N=4, L=10, H=1)],
)
def test_invalid_ring(kwargs):
    with pytest.raises(ValueError):
        Ring(**kwargs)


def test_ratio_override_and_gamma_range():
    eq = equivalent_2d(Ring(**SYNTHETIC, ratio_override=2.0))
    assert eq.ratio_overridden and eq.gamma_eff(0.1) == pytest.approx(0.2)
    with pytest.raises(ValueError):
        eq.gamma_eff(1.5)


def test_cli(tmp_path):
    out = tmp_path / "out.csv"
    main([str(EXAMPLES / "geometry_example.csv"), "--reactions",
          str(EXAMPLES / "reactions_example.csv"), "-o", str(out)])
    rows = list(csv.reader(out.open()))
    assert len(rows) == 7 and rows[0][-2].startswith("gamma_eff[R1")
