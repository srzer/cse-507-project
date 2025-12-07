from dataclasses import dataclass
from math import prod
from typing import Iterable

from dreal import And, Box, Formula, Variable
from typing_extensions import Self

from interval import Bounds

from .point_type import Point


# immutable
@dataclass(frozen=True)
class BoxN:
    min: Point
    max: Point

    def __post_init__(self: Self):
        assert len(self.min) == len(self.max), "Box dimension mismatch."

    @classmethod
    def from_lists(cls, l1: Iterable[float], l2: Iterable[float]) -> Self:
        return cls(Point(l1), Point(l2))

    def __str__(self: Self) -> str:
        return f"({self.min}, {self.max})"

    @property
    def dim(self: Self) -> int:
        return len(self.min)

    @property
    def sides(self: Self) -> list[Bounds]:
        return [Bounds(lo, hi) for lo, hi in zip(self.min, self.max)]

    # def sides(self: Self) -> list[Tuple[float, float]]: return list(zip(self.min, self.max))

    @property
    # center point of box
    def center(self: Self) -> Point:
        return Point(tuple((hi + lo) / 2 for lo, hi in self.sides))

    @property
    # lengths of all sides
    def lengths(self: Self) -> list[float]:
        return [hi - lo for lo, hi in self.sides]

    # center of side
    def _side_center(self: Self, i: int):
        assert 0 <= i <= self.dim, f"Side index {i} our of range {self.dim}."
        return (self.max[i] + self.min[i]) / 2

    @property
    # length of longest side
    def max_side_length(self: Self) -> float:
        return max(self.lengths)

    @property
    # midpoint of longest side
    def max_side_center(self: Self) -> float:
        return self._side_center(self._max_side_idx())

    @property
    # calculate box volume
    def volume(self: Self) -> float:
        return prod(self.lengths)

    # return index of dimension with longest side
    def _max_side_idx(self: Self) -> int:
        return self.lengths.index(max(self.lengths))

    # check if box containts a point
    def contains(self: Self, p: Point) -> bool:
        assert len(p) == self.dim, (
            f"Point dimension {len(p)} does not match Box dimension {self.dim}"
        )
        return all(lo <= pt <= hi for pt, lo, hi in zip(p, self.min, self.max))

    # build dReal constraints for box
    def build_constraints(self: Self, vars: list[Variable]) -> Formula:
        if len(vars) != self.dim:
            raise ValueError(f"Expected {self.dim} variables, got {len(vars)}")
        # this is some polymorphic insanity...
        constraints = [
            constr
            for var, lo, hi in zip(vars, self.min, self.max)
            for constr in (lo <= var, var <= hi)
        ]
        return And(*constraints) if constraints else Formula.TRUE()


# convert from dReal Box type
def from_box_model(model: Box) -> BoxN:
    ivs = model.values()
    return BoxN.from_lists([iv.lb() for iv in ivs], [iv.ub() for iv in ivs])


# build initial axis-aligned n dimensional box
def build_basic_box(min, max: float, dim: int) -> BoxN:
    return BoxN(Point([min] * dim), Point([max] * dim))
