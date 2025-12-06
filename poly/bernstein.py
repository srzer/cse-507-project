from itertools import product
from dataclasses import dataclass, field
from collections import UserDict
from math import comb
from typing import Dict, Tuple, NewType, Iterator
from typing_extensions import Self

from interval import Bounds
from box import BoxN
from poly.type import Polynomial, Rational, to_term

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


ZERO_TOLERANCE = 1e-15  # tolerance for filtering zero coefficients
DENOM_MIN_THRESHOLD = 1e-8  # min demoninator to prevent div-by-zero


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


# old bernstein code
# def bernstein_coefficients_nd(poly: Polynomial, box: BoxN, degrees=None):
#     """
#     Correct multivariate Bernstein coefficient computation on a box.

#     poly:   Polynomial (list of (coeff, exponents) terms) in power basis.
#     box:    list/tuple of (a_j, b_j) intervals.
#     degrees:
#         Optional list of degrees per variable. If None, the maximal exponent
#         in poly for each variable is used.

#     Returns:
#         (B, degrees), where B maps multi-indices (k1,...,kn) to Bernstein
#         coefficients on the given box, and degrees is the degree per variable.
#     """
#     n = len(box)

#     # 1) Determine degrees per variable
#     if degrees is None:
#         degrees = [0] * n
#         for coeff, mon in poly:
#             for j in range(n):
#                 degrees[j] = max(degrees[j], mon[j])

#     # 2) Build a sparse tensor of power-basis coefficients in t after
#     #    the change of variables x_j = a_j + (b_j - a_j) * t_j.
#     #    P[k1,...,kn] is the coefficient of t_1^k1 * ... * t_n^kn.
#     P = {}

#     for coeff, mon in poly:
#         # mon is the exponents in x
#         # For each dimension j, expand x_j^{mon[j]} into t_j^{k_j}:
#         #   x_j^{m_j} = sum_{k_j=0}^{m_j} comb(m_j, k_j) * a_j^{m_j-k_j} * (b_j-a_j)^{k_j} * t_j^{k_j}
#         ranges = [range(mj + 1) for mj in mon]
#         for ks in product(*ranges):
#             term_coeff = coeff
#             for j, kj in enumerate(ks):
#                 mj = mon[j]
#                 a_j, b_j = box[j]
#                 term_coeff *= comb(mj, kj) * (a_j ** (mj - kj)) * ((b_j - a_j) ** kj)
#             P[ks] = P.get(ks, 0.0) + term_coeff

#     # 3) Convert the power-basis tensor P into a Bernstein-basis tensor B
#     #    along each axis using the 1D transformation:
#     #
#     #   For a univariate polynomial of degree d,
#     #       p(t) = sum_{i=0}^d a_i t^i
#     #   the Bernstein coefficients {b_k} of degree d on [0,1] satisfy:
#     #       b_k = sum_{i=0}^k comb(k, i) / comb(d, i) * a_i.
#     #
#     B = dict(P)

#     for axis in range(n):
#         d = degrees[axis]
#         newB = {}
#         # Fix all other indices, and transform the coefficients along this axis.
#         other_axes = [i for i in range(n) if i != axis]
#         other_ranges = [range(degrees[i] + 1) for i in other_axes]

#         for other_idx in product(*other_ranges):
#             # Collect a_i for i = 0..d along this axis
#             a_vec = []
#             for i in range(d + 1):
#                 idx = []
#                 oi = 0
#                 for ax in range(n):
#                     if ax == axis:
#                         idx.append(i)
#                     else:
#                         idx.append(other_idx[oi])
#                         oi += 1
#                 idx = tuple(idx)
#                 a_vec.append(B.get(idx, 0.0))

#             # Power -> Bernstein conversion
#             b_vec = [0.0] * (d + 1)
#             for k in range(d + 1):
#                 s = 0.0
#                 for i in range(k + 1):
#                     if comb(d, i) == 0:
#                         continue
#                     s += comb(k, i) / comb(d, i) * a_vec[i]
#                 b_vec[k] = s

#             # Store back into newB
#             for k, val in enumerate(b_vec):
#                 idx = []
#                 oi = 0
#                 for ax in range(n):
#                     if ax == axis:
#                         idx.append(k)
#                     else:
#                         idx.append(other_idx[oi])
#                         oi += 1
#                 newB[tuple(idx)] = val

#         B = newB

#     return B, degrees


# def bernstein_bounds(B: Dict[Tuple[int, ...], float]) -> Tuple[float, float]:
#     """Return min and max of a Bernstein coefficient tensor B (dict)."""
#     vals = list(B.values())
#     return min(vals), max(vals)


# def bernstein_bounds_on_box(
#     poly_num: Polynomial, poly_den: Polynomial, box: BoxND
# ) -> Tuple[float, float]:
#     """
#     Compute a crude lower/upper bound for f = N/D on a BoxND using Bernstein.

#     This uses your existing bernstein_coefficients_nd and bernstein_bounds
#     utilities. It ignores additional constraints; constraints are still
#     handled by dReal in the main algorithm.

#     Args:
#         poly_num, poly_den: Polynomial (list of (coeff, exponents) terms)
#         box: BoxND (lows, highs)

#     Returns:
#         (lb, ub): lower and upper bounds from Bernstein coefficients.
#     """
#     # Convert BoxND to a list of [a_i, b_i] intervals
#     intervals = list(zip(box.lows, box.highs))
#     n = len(intervals)

#     # Determine global degrees for each variable
#     degrees = [0] * n

#     def update_degrees(poly: Polynomial):
#         for coeff, mon in poly:
#             for j in range(n):
#                 degrees[j] = max(degrees[j], mon[j])

#     update_degrees(poly_num)
#     update_degrees(poly_den)

#     # Compute Bernstein coefficients of numerator and denominator
#     B_num, _ = bernstein_coefficients_nd(poly_num, intervals, degrees=degrees)
#     B_den, _ = bernstein_coefficients_nd(poly_den, intervals, degrees=degrees)

#     # Get min/max of Bernstein coefficients
#     N_min, N_max = bernstein_bounds(B_num)
#     D_min, D_max = bernstein_bounds(B_den)

#     # Protect the denominator from being too small
#     D_min = max(D_min, 1e-8)
#     D_max = max(D_max, 1e-8)

#     lb = N_min / D_max
#     ub = N_max / D_min
#     # print(intervals, N_min, N_max, D_min, D_max, lb)
#     return lb, ub


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
def bernstein_bounds(f: Rational, box: BoxN) -> Bounds:
    num_min, num_max = bernstein_coefficients(f.num, box.sides).bounds
    den_min, den_max = bernstein_coefficients(f.den, box.sides).bounds

    # protect denominator
    den_min, den_max = (
        max(den_min, DENOM_MIN_THRESHOLD),
        max(den_max, DENOM_MIN_THRESHOLD),
    )

    return Bounds(num_min / den_max, num_max / den_min)
