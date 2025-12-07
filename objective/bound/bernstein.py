from collections import UserDict
from dataclasses import dataclass, field
from itertools import product
from math import comb
from typing import Dict, Iterator, NewType, Tuple

from typing_extensions import Self

from box import BoxN
from interval import Bounds

from objective.polynomial_type import Polynomial, to_term
from objective.rational_type import Rational

from .type import ObjectiveBounds, ZERO_TOLERANCE, DENOM_MIN_THRESHOLD


class BernsteinBounds(ObjectiveBounds):
    def _run(self: Self, obj: Rational, box: BoxN) -> Bounds:
        return _bernstein_bounds(obj, box)


# {[exponents]: coeff}
MonomialDict = NewType("MonomialDict", Dict[Tuple[int, ...], float])


@dataclass
# multi-index to coeff mapping for Bernstein basis
class BernsteinCoeffs(UserDict[Tuple[int, ...], float]):
    data: dict[Tuple[int, ...], float] = field(default_factory=dict)

    def __getitem__(self, key: Tuple[int, ...]) -> float:
        return self.data[key]

    def __setitem__(self, key: Tuple[int, ...], value: float) -> None:
        self.data[key] = value

    def items(self):
        return self.data.items()

    def values(self):
        return self.data.values()

    def __iter__(self) -> Iterator[Tuple[int, ...]]:
        return iter(self.data)

    @property
    # return min and max of coefficients
    def bounds(self: Self) -> Bounds:
        if not self.data:
            return Bounds(0.0, 0.0)  # NOTE: default singleton zero interval
        vals = list(self.data.values())
        return Bounds(min(vals), max(vals))


# convert to dict, combining like terms
def to_monomial_dict(poly: Polynomial) -> MonomialDict:
    result = {}
    for coeff, exps in poly:
        result[exps] = result.get(exps, 0.0) + coeff
    return MonomialDict(result)


# convert from monomial dict, filtering zero coefficients
def from_monomial_dict(d: MonomialDict) -> Polynomial:
    return Polynomial(
        [
            to_term(coeff, list(exps))
            for exps, coeff in d.items()
            if abs(coeff) > ZERO_TOLERANCE
        ]
    )


# conver 1-dim power basis coeffs to Bernstein basis
def power_to_bernstein_1d(coeffs: list[float]) -> list[float]:
    d = len(coeffs) - 1
    bernstein = [0.0] * (d + 1)

    for k in range(d + 1):
        binom_d_k = comb(d, k)
        if binom_d_k == 0:
            continue
        bernstein[k] = sum(
            comb(k, i) / comb(d, i) * coeffs[i] for i in range(k + 1) if comb(d, i) != 0
        )

    return bernstein


def transform_to_unit_box(
    poly: MonomialDict, intervals: list[Bounds]
) -> BernsteinCoeffs:
    """
    Transform polynomial from [a,b]^n to [0,1]^n via x_j = a_j + (b_j - a_j)*t_j.
    Returns coefficients in the t variables.
    """
    result = {}
    for exps, coeff in poly.items():
        # expand each x_j^m_j into sum over t_j^k_j
        for ks in product(*(range(m + 1) for m in exps)):
            term_coeff = coeff
            for j, (k, m) in enumerate(zip(ks, exps)):
                a, b = intervals[j]
                term_coeff *= comb(m, k) * (a ** (m - k)) * ((b - a) ** k)
            result[ks] = result.get(ks, 0.0) + term_coeff

    return BernsteinCoeffs(result)


# apply power to Bernstein transformation along a single axis
def apply_bernstein_transform_axis(
    coeffs: BernsteinCoeffs, axis: int, degrees: list[int]
) -> BernsteinCoeffs:
    result = BernsteinCoeffs({})

    # get indices for other axes
    other_ranges = [range(degrees[i] + 1) for i in range(len(degrees)) if i != axis]

    for other_idx in product(*other_ranges):
        # extract 1D slice along this axis
        slice_coeffs = []
        for i in range(degrees[axis] + 1):
            idx = list(other_idx)
            idx.insert(axis, i)
            slice_coeffs.append(coeffs.get(tuple(idx), 0.0))

        # transform this slice
        bernstein_slice = power_to_bernstein_1d(slice_coeffs)

        # store back
        for k, val in enumerate(bernstein_slice):
            idx = list(other_idx)
            idx.insert(axis, k)
            result[tuple(idx)] = val

    return result


# compute multivariable Bernstein coefficients on a box
def bernstein_coefficients(
    poly: Polynomial,
    intervals: list[Bounds],
) -> BernsteinCoeffs:
    """
    Args:
        poly: Polynomial in power basis
        intervals: Box as list of (low, high) per dimension
    Returns:
        (coeffs, degrees) where coeffs maps multi-indices to Bernstein coefficients
    """
    degrees = poly.max_var_degrees()

    # convert to dict and transform to unit box
    poly_dict = to_monomial_dict(poly)
    power_coeffs = transform_to_unit_box(poly_dict, intervals)

    # apply Bernstein transformation along each axis
    bernstein_coeffs = power_coeffs

    for axis in range(len(intervals)):
        bernstein_coeffs = apply_bernstein_transform_axis(
            bernstein_coeffs, axis, degrees
        )

    return bernstein_coeffs


# compute bounds for ration function f on a box using Bernstein basis.
def _bernstein_bounds(f: Rational, box: BoxN) -> Bounds:
    num_min, num_max = bernstein_coefficients(f.num, box.sides).bounds
    den_min, den_max = bernstein_coefficients(f.den, box.sides).bounds

    # protect denominator
    den_min, den_max = (
        max(den_min, DENOM_MIN_THRESHOLD),
        max(den_max, DENOM_MIN_THRESHOLD),
    )

    return Bounds(num_min / den_max, num_max / den_min)
