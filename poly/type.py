import math
from dataclasses import dataclass
from typing import Any, List, NewType, Tuple

from dreal import Expression, Variable

from box import Point

"""
     exponent list length = # of vars
    [coeff: [exp1, exp2, ..., expN]]
"""

Term = NewType("Term", Tuple[float, Tuple[int, ...]])  # (coeff, [exponents])
Polynomial = NewType("Polynomial", List[Term])

# MonomialDict = Dict[Tuple[int, ...], float] # (coeff, [exponents])


@dataclass(frozen=True)
class Rational:
    num: Polynomial
    den: Polynomial


def to_term(coeff: float, exps: list[int]) -> Term:
    return Term((coeff, tuple(exps)))


# def to_monomial_dict(poly: Polymonial) -> MonomialDict


# kinda cool that this works so polymorphically
def _eval(poly: Polynomial, xs: list[Any]) -> Any:
    return sum(
        coeff * math.prod(x_j**e for x_j, e in zip(xs, mon) if e != 0)
        for coeff, mon in poly
    )


# evaluate polynomial numerically at point p
def eval(f: Rational, p: Point) -> float:
    return float(_eval(f.num, list(p)) / _eval(f.den, list(p)))


# evaluate polynomial symbolically with variable list xs
def eval_symbolic(f: Rational, xs: list[Variable]) -> Expression:
    return _eval(f.num, xs) / _eval(f.den, xs)
