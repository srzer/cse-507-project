from typing_extensions import Self

from box import BoxN
from interval import AffineInterval, Bounds

from objective.polynomial_type import _eval_polynomial_affine_bounds
from objective.rational_type import Rational

from .type import ObjectiveBounds, ZERO_TOLERANCE


class AffineBounds(ObjectiveBounds):
    def _run(self: Self, obj: Rational, box: BoxN) -> Bounds:
        return _affine_bounds(obj, box)


def _affine_bounds(
    obj: Rational,
    box: BoxN,
) -> Bounds:
    """
    Compute a cheap affine/interval-based lower/upper bound for

        f(x) = N(x) / D(x)

    over a BoxND.

    This uses very lightweight affine forms: each variable is represented as
    a midpoint plus a single noise term whose range matches the box interval.
    We then evaluate numerator and denominator bounds using eval_poly_aa.

    Args:
        poly_num, poly_den:
            monomial dicts for numerator and denominator, respectively.
        box:
            A BoxND-like object with attributes:
                - lows:  list/sequence of lower bounds
                - highs: list/sequence of upper bounds
        safe_eps:
            Used to clip the denominator away from zero (and negative values).
    Returns:
        (lb, ub):
            lower and upper bounds on f over this box.
    """
    # build one affine variable per dimension, covering the entire box interval.
    affine_vars = [
        AffineInterval(mid, [diam / 2]) for mid, diam in zip(box.center, box.lengths)
    ]

    # bounds for numerator and denominator
    num_min, num_max = _eval_polynomial_affine_bounds(obj.num, affine_vars)
    den_min, den_max = _eval_polynomial_affine_bounds(obj.den, affine_vars)

    # protect denominator from being too small or negative (assuming den > 0 on feasible set)
    den_min = max(den_min, ZERO_TOLERANCE)
    den_max = max(den_max, ZERO_TOLERANCE)

    # rational bounds: for den > 0,
    #   f_min >= num_min / den_max
    #   f_max <= num_max / den_min
    return Bounds(num_min / den_max, num_max / den_min)
