import math
from dataclasses import dataclass, replace
from typing import Tuple
from typing_extensions import Self

from dreal import Box


# Point = NewType("Point", Tuple[float, ...])
@dataclass(frozen=True)
class Point:
    coords: Tuple[float, ...]

    def with_value(self, i: int, v: float) -> Self:
        assert i < len(self), f"Point index {i} out of range {len(self)}."
        new = self.coords[:i] + (v,) + self.coords[i + 1 :]
        return replace(self, coords=Point(new))

    def __getitem__(self, i: int) -> float:
        return self.coords[i]

    def __iter__(self):
        return iter(self.coords)

    def __len__(self):
        return len(self.coords)


# immutable
@dataclass(frozen=True)
class BoxN:
    min: Point
    max: Point

    # this triggers at some point...?
    def __post_init__(self: Self):
        assert len(self.min) == len(self.max), "Box dimension mismatch."

    @property
    def dim(self: Self):
        return len(self.min)

    @property
    # center point of box
    def center(self: Self) -> Point:
        return Point(tuple((hi + lo) / 2 for lo, hi in zip(self.min, self.max)))

    @property
    # lengths of all sides
    def lengths(self: Self) -> list[float]:
        return [hi - lo for lo, hi in zip(self.min, self.max)]

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
        sides = [hi - lo for lo, hi in zip(self.min, self.max)]
        return math.prod(sides)

    # return index of dimension with longest side
    def _max_side_idx(self: Self) -> int:
        sides = [hi - lo for lo, hi in zip(self.min, self.max)]
        return sides.index(max(sides))

    # check if box containts a point
    def contains(self: Self, point: Point) -> bool:
        if len(point) != self.dim:
            raise ValueError(
                f"Point dimension {len(point)} does not match Box dimension {self.dim}"
            )
        return all(lo <= pt <= hi for pt, lo, hi in zip(point, self.min, self.max))


def from_box_model(model: Box) -> BoxN:
    return BoxN(
        Point(tuple([p.lb() for p in model])), Point(tuple([p.ub() for p in model]))
    )
