import math
from dataclasses import dataclass
from typing import Any, NewType, Tuple

from typing_extensions import Self

from box import Point
from interval import AffineInterval, Bounds

# exponent list length = # of vars
Term = NewType("Term", Tuple[float, Tuple[int, ...]])  # (coeff, [exponents])


@dataclass(frozen=True)
# polynomial in power basis form
# all terms must have exponent tuples of same length
class Polynomial:
    terms: list[Term]

    # validate all exponent tuples have consistent dimension.
    def __post_init__(self):
        if not self.terms:
            return

        # use expected dimension from first term
        n = len(self.terms[0][1])

        # validate all terms match this dimension
        for coeff, exps in self.terms:
            if len(exps) != n:
                raise ValueError(
                    f"Inconsistent exponent dimensions in polynomial: "
                    f"expected {n} variables, got {len(exps)} in term ({coeff}, {exps})"
                )

    @property
    # number of vars in polynomial
    def n_vars(self) -> int:
        return len(self.terms[0][1]) if self.terms else 0

    # get max degree for each vav
    def max_var_degrees(self: Self) -> list[int]:
        if not self.terms:
            return []

        degrees = [0] * self.n_vars
        for _, exps in self.terms:
            for i, exp in enumerate(exps):
                degrees[i] = max(degrees[i], exp)
        return degrees

    def __iter__(self):
        return iter(self.terms)


# constant unit polynomial = 1
def unit(dim: int) -> Polynomial:
    return Polynomial([to_term(1, [0] * dim)])


# create a term from coeff and exponent list
def to_term(coeff: float, exps: list[int]) -> Term:
    return Term((coeff, tuple(exps)))


# create a polynomial with expontent length validation
def to_polynomial(terms: list[Term]) -> Polynomial:
    return Polynomial(terms)


# kinda cool that this works so polymorphically
def _eval_polynomial(poly: Polynomial, xs: list[Any]) -> Any:
    return sum(
        coeff * math.prod(x_j**e for x_j, e in zip(xs, mon) if e != 0)
        for coeff, mon in poly
    )


def _diff_poly(f: Polynomial, p: Point) -> list[float]:
    n = len(p)
    grad = [0.0] * n

    for coeff, mon in f:
        # for each variable k, compute ∂/∂x_k of this term
        for k in range(n):
            e_k = mon[k]
            if e_k == 0:
                continue
            # power rule: ∂/∂x_k (c * x_k^e_k * ...) = c * e_k * x_k^(e_k-1) * ...
            term_grad = coeff * e_k * (p[k] ** (e_k - 1))
            for j in range(n):
                if j == k:
                    continue
                term_grad *= p[j] ** mon[j]
            grad[k] += term_grad

    return grad


def _eval_polynomial_affine_bounds(
    poly: Polynomial, affine_vars: list[AffineInterval]
) -> Bounds:
    """
    Evaluate bounds for a polynomial over a box using interval arithmetic.

    Args:
        poly: Polynomial in power basis form
        affine_vars: List of AffineInterval objects, one per variable

    Returns:
        Bounds object with min and max values
    """
    min_val = 0.0
    max_val = 0.0

    for coeff, exponents in poly.terms:
        # start with coefficient bounds
        term_min = term_max = coeff

        # for each variable, multiply by x_j^e bounds
        for j, e in enumerate(exponents):
            if e == 0:
                continue  # x^0 = 1, no effect

            # get variable bounds
            bounds = affine_vars[j].to_bounds()
            a_j, b_j = bounds.lo, bounds.hi

            # compute x_j^e bounds on [a_j, b_j]
            if e < 0:
                raise ValueError(f"Negative exponent {e} not supported for polynomials")

            if a_j >= 0:
                # monotonic on positive interval
                p_min, p_max = a_j**e, b_j**e
            elif b_j <= 0:
                # monotonic on negative interval
                p_min, p_max = b_j**e, a_j**e
            else:
                # interval crosses zero
                p_min = 0.0
                p_max = max(abs(a_j), abs(b_j)) ** e

            # update term bounds
            new_bounds = [
                term_min * p_min,
                term_min * p_max,
                term_max * p_min,
                term_max * p_max,
            ]
            term_min, term_max = min(new_bounds), max(new_bounds)

        min_val += term_min
        max_val += term_max

    return Bounds(min_val, max_val)
