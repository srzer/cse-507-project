from dataclasses import dataclass

from dreal import Expression, Variable
from typing_extensions import Self

from box import Point
from objective.polynomial_type import Polynomial, _eval_polynomial, _diff_poly


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


# evaluate rational function numerically at point p
def eval_rational(f: Rational, p: Point) -> float:
    return float(_eval_polynomial(f.num, list(p)) / _eval_polynomial(f.den, list(p)))


# evaluate rational function symbolically with variable list xs
def eval_symbolic(f: Rational, xs: list[Variable]) -> Expression:
    return _eval_polynomial(f.num, xs) / _eval_polynomial(f.den, xs)


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
