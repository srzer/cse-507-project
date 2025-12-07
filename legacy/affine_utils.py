# affine_utils.py

"""
Lightweight affine / interval helpers for bounding a rational function

    f(x) = N(x) / D(x)

on an axis-aligned box.

This module is designed to mirror the structure of the Bernstein utilities:
- we provide a cheap bounding routine affine_bounds_on_box(poly_num, poly_den, box)
- polynomials are given as monomial dicts: {(e1,...,en): coeff}
- 'box' is assumed to be a BoxND-like object with .lows and .highs
"""

from typing import Dict, List, Tuple


class AffineForm:
    """
    Very lightweight affine form:

        x ≈ mid + coeffs[0] * ε,   with ε in [-1, 1].

    We only use it to get a simple range via:

        [mid - |coeffs[0]|, mid + |coeffs[0]|]

    No full affine propagation or correlation tracking is done.
    """

    def __init__(self, mid: float, coeffs: List[float] = None):
        self.mid = mid
        self.coeffs = coeffs if coeffs is not None else []

    def range(self) -> Tuple[float, float]:
        """
        Return the interval [min, max] represented by this affine form.
        """
        r = sum(abs(c) for c in self.coeffs)
        return (self.mid - r, self.mid + r)


def eval_poly_aa(poly: Dict[Tuple[int, ...], float],
                 affine_vars: List[AffineForm]) -> Tuple[float, float]:
    """
    Evaluate rough bounds [min_val, max_val] for a polynomial over a box
    using only the ranges of each variable (effectively interval arithmetic).

    poly:
        dict {(e1,...,en): coeff}  in power basis
    affine_vars:
        list[AffineForm], one per variable

    Returns:
        (min_val, max_val)
    """
    n = len(affine_vars)
    min_val = 0.0
    max_val = 0.0

    for mon, coeff in poly.items():
        # Start with the coefficient
        term_min = coeff
        term_max = coeff

        for j in range(n):
            a_j, b_j = affine_vars[j].range()
            e = mon[j]

            # Evaluate x_j^e on [a_j, b_j]
            if e == 0:
                # x^0 = 1
                p_min = p_max = 1.0
            elif e > 0:
                # If interval does not cross 0, x^e is monotone.
                # If it crosses 0, min is 0 and max is max(|a|,|b|)^e.
                if a_j >= 0:
                    p_min = a_j ** e
                    p_max = b_j ** e
                elif b_j <= 0:
                    p_min = b_j ** e
                    p_max = a_j ** e
                else:
                    # crosses zero
                    p_min = 0.0
                    p_max = max(abs(a_j), abs(b_j)) ** e
            else:
                # Negative exponents are not expected for polynomials.
                raise ValueError("Negative exponent in polynomial monomial")

            # Multiply current term bounds by [p_min, p_max]
            candidates = [
                term_min * p_min,
                term_min * p_max,
                term_max * p_min,
                term_max * p_max,
            ]
            term_min = min(candidates)
            term_max = max(candidates)

        min_val += term_min
        max_val += term_max

    return (min_val, max_val)


def affine_bounds_on_box(poly_num: Dict[Tuple[int, ...], float],
                         poly_den: Dict[Tuple[int, ...], float],
                         box,
                         safe_eps: float = 1e-8) -> Tuple[float, float]:
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
    # Build one affine variable per dimension, covering the entire box interval.
    affine_vars: List[AffineForm] = []
    for lo, hi in zip(box.lows, box.highs):
        mid = 0.5 * (lo + hi)
        rad = 0.5 * (hi - lo)
        affine_vars.append(AffineForm(mid, [rad]))

    # Bounds for numerator and denominator
    num_min, num_max = eval_poly_aa(poly_num, affine_vars)
    den_min, den_max = eval_poly_aa(poly_den, affine_vars)

    # Protect denominator from being too small or negative (assuming D > 0 on feasible set)
    den_min = max(den_min, safe_eps)
    den_max = max(den_max, safe_eps)

    # Rational bounds: for D > 0,
    #   f_min >= num_min / den_max
    #   f_max <= num_max / den_min
    lb = num_min / den_max
    ub = num_max / den_min

    return lb, ub
