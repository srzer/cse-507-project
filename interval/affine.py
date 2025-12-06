from dataclasses import dataclass
from typing_extensions import Self

from .type import Bounds

"""
Basic affine / interval helpers for bounding a rational function

    f(x) = N(x) / D(x)

on an axis-aligned box.

This module is designed to mirror the structure of the Bernstein utilities:
- we provide a cheap bounding routine affine_bounds_on_box(poly_num, poly_den, box)
- polynomials are given as monomial dicts: {(e1,...,en): coeff}
- 'box' is assumed to be a BoxND-like object with .lows and .highs
"""


@dataclass(frozen=True)
class AffineInterval:
    """
    simple 'affine' form:
      x ≈ mid + sum_i coeffs[i] * ε_i,  with ε_i ∈ [-1,1].
    We are not doing full affine propagation here; we only use
    mid/coeffs to get a cheap range and rescale for sub-boxes.
    We only use it to get a simple range via:
        [mid - |coeffs[0]|, mid + |coeffs[0]|]
    No full affine propagation or correlation tracking is done.
    """

    mid: float
    coeffs: list[float]

    # return interval [min, max] represented by affine form
    def to_bounds(self: Self) -> Bounds:
        r = sum(abs(c) for c in self.coeffs)
        return Bounds(self.mid - r, self.mid + r)
