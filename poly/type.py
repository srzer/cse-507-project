import math
from dataclasses import dataclass
from typing import Any, List, NewType, Tuple

from dreal import Expression, Variable
from typing_extensions import Self

from box import Point

# exponent list length = # of vars
Term = NewType("Term", Tuple[float, Tuple[int, ...]])  # (coeff, [exponents])


@dataclass(frozen=True)
# polynomial in power basis form
# all terms must have exponent tuples of same length
class Polynomial:
    terms: List[Term]

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


@dataclass(frozen=True)
class Rational:
    num: Polynomial
    den: Polynomial

    # validate num and den have same # of vars
    def __post_init__(self: Self):
        if self.num.n_vars != self.den.n_vars:
            raise ValueError(
                f"Numerator and denominator must have same dimensions: "
                f"num has {self.num.n_vars} vars, den has {self.den.n_vars} vars"
            )

    @property
    # number of variables in rational function
    def n_vars(self) -> int:
        return self.num.n_vars


# create a term from coeff and exponent list
def to_term(coeff: float, exps: list[int]) -> Term:
    return Term((coeff, tuple(exps)))


# create a polynomial with expontent length validation
def to_polynomial(terms: List[Term]) -> Polynomial:
    return Polynomial(terms)


# kinda cool that this works so polymorphically
def _eval_polynomial(poly: Polynomial, xs: list[Any]) -> Any:
    return sum(
        coeff * math.prod(x_j**e for x_j, e in zip(xs, mon) if e != 0)
        for coeff, mon in poly
    )


# evaluate rational function numerically at point p
def eval_rational(f: Rational, p: Point) -> float:
    return float(_eval_polynomial(f.num, list(p)) / _eval_polynomial(f.den, list(p)))


# evaluate rational function symbolically with variable list xs
def eval_symbolic(f: Rational, xs: list[Variable]) -> Expression:
    return _eval_polynomial(f.num, xs) / _eval_polynomial(f.den, xs)


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


# evaluate the gradient vector of the given function f at point p
def eval_gradient(f: Rational, p: Point) -> list[float]:
    n = len(p)

    num = _eval_polynomial(f.num, list(p))
    den = _eval_polynomial(f.den, list(p))

    grad_num = _diff_poly(f.num, p)
    grad_den = _diff_poly(f.den, p)

    # prevent division by zero
    den = max(abs(den), 1e-12)

    # quotient rule: ∇(u/v) = (∇u * v - u * ∇v) / v²
    grad_f = []
    for k in range(n):
        grad_f_k = (grad_num[k] * den - num * grad_den[k]) / (den**2)
        grad_f.append(grad_f_k)

    return grad_f
