"""3D-to-2D equivalence of multi-hole showerhead rings (Yoon & Lee, Nature Chemical Engineering, submitted).

Flow-resistance equivalence (article Eq. (3)): the equivalent 2D annular
half-width R' is the root of

    N r^4 = R_o^4 [ (1 - k^4) - (1 - k^2)^2 / ln(1/k) ]

with R_o = L + R', R_i = L - R', k = R_i / R_o, i.e. the summed Poiseuille
conductance of the N holes equals the conductance of the annulus.

Reaction equivalence (article Eq. (4)):

    gamma_eff = gamma * A_wall,3D / A_wall,2D

applied to the hole side-walls only (chamber walls and heater keep gamma).

All lengths are in mm. The numerics (series fallback, 200-step bisection)
match the web calculator at https://junghwanyoon.dev/tools/showerhead-calculator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Below this 1 - k the closed-form bracket loses precision to cancellation.
SERIES_EPS = 5e-4
BISECTION_STEPS = 200


def annulus_bracket(k: float) -> float:
    """(1 - k^4) - (1 - k^2)^2 / ln(1/k), with a series fallback near k -> 1."""
    eps = 1 - k
    if eps < SERIES_EPS:
        # Taylor expansion around k = 1: (4/3) eps^3 (1 - eps/2) + O(eps^5)
        return (4 / 3) * eps**3 * (1 - eps / 2)
    return 1 - k**4 - (1 - k * k) ** 2 / math.log(1 / k)


def mismatch(r_prime: float, L: float, N: int, r: float) -> float:
    """Conductance mismatch g(R'); monotonically increasing in R'."""
    ro = L + r_prime
    k = (L - r_prime) / ro
    return ro**4 * annulus_bracket(k) - N * r**4


def solve_r_prime(N: int, r: float, L: float) -> float:
    """Solve Eq. (3) for R' by bisection on (0, L).

    Raises ValueError when N r^4 exceeds the largest annulus conductance
    available at this L (no equivalent annulus exists).
    """
    lo = 0.0
    hi = L * (1 - 1e-9)
    if mismatch(hi, L, N, r) < 0:
        raise ValueError(
            "no equivalent annulus exists (N*r^4 exceeds the maximum annulus "
            "conductance for this L) - reduce N or r, or increase L"
        )
    for _ in range(BISECTION_STEPS):
        mid = (lo + hi) / 2
        if mismatch(mid, L, N, r) >= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


@dataclass(frozen=True)
class Ring:
    """One ring of N holes of radius r at radial position L, plate thickness H [mm]."""

    r: float
    N: int
    L: float
    H: float
    label: str = ""
    ratio_override: float | None = None

    def __post_init__(self) -> None:
        if not (self.r > 0 and self.L > 0 and self.H > 0):
            raise ValueError(f"{self.label or 'ring'}: r, L and H must be positive")
        if not (isinstance(self.N, int) and self.N > 0):
            raise ValueError(f"{self.label or 'ring'}: N must be a positive integer")
        if self.r >= self.L:
            raise ValueError(f"{self.label or 'ring'}: r must be smaller than L")
        if self.ratio_override is not None and not self.ratio_override > 0:
            raise ValueError(f"{self.label or 'ring'}: ratio_override must be positive")


@dataclass(frozen=True)
class Equivalent2D:
    """Equivalent annulus of a ring and its diagnostics."""

    r_prime: float
    ro: float
    ri: float
    k: float
    open_area_3d: float
    open_area_2d: float
    wall_area_3d: float
    wall_area_2d: float
    area_ratio: float
    ratio_overridden: bool
    residual_rel: float
    """Relative residual of Eq. (3) at the returned R' (root-finding accuracy)."""
    r_prime_error_est: float
    """Residual converted to an R' error estimate via dg/dR' [mm]."""

    def gamma_eff(self, gamma: float) -> float:
        """Area-corrected reaction probability for the hole side-walls (Eq. (4))."""
        if not 0 < gamma <= 1:
            raise ValueError("gamma must be in (0, 1]")
        return gamma * self.area_ratio


def equivalent_2d(ring: Ring) -> Equivalent2D:
    """Equivalent 2D annulus of a hole ring (Eq. (3)) and its wall-area ratio."""
    N, r, L, H = ring.N, ring.r, ring.L, ring.H
    r_prime = solve_r_prime(N, r, L)
    ro = L + r_prime
    ri = L - r_prime
    wall_area_3d = N * 2 * math.pi * r * H
    wall_area_2d = 2 * math.pi * H * (ro + ri)

    # Solver diagnostics: residual_rel = g(R') / (N r^4);
    # r_prime_error_est = |g(R')| / |dg/dR'| (Newton step, first-order bound).
    target = N * r**4
    residual = mismatch(r_prime, L, N, r)
    delta = min(L * 1e-7, r_prime / 2, (L - r_prime) / 2)
    dg = (
        (mismatch(r_prime + delta, L, N, r) - mismatch(r_prime - delta, L, N, r)) / (2 * delta)
        if delta > 0
        else 0.0
    )

    return Equivalent2D(
        r_prime=r_prime,
        ro=ro,
        ri=ri,
        k=ri / ro,
        open_area_3d=N * math.pi * r**2,
        open_area_2d=math.pi * (ro**2 - ri**2),
        wall_area_3d=wall_area_3d,
        wall_area_2d=wall_area_2d,
        area_ratio=ring.ratio_override if ring.ratio_override is not None else wall_area_3d / wall_area_2d,
        ratio_overridden=ring.ratio_override is not None,
        residual_rel=residual / target,
        r_prime_error_est=abs(residual / dg) if dg != 0 else 0.0,
    )
