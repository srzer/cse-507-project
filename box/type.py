import math
from dataclasses import dataclass
from typing import NewType, Tuple, Self

from dreal import Box, Interval, Variable

Point = NewType("Point", Tuple[float, ...])


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
    # get center point of box
    def center(self: Self) -> Point:
        return Point(tuple((hi + lo) / 2 for lo, hi in zip(self.min, self.max)))

    @property
    # return length of longest side
    def max_side(self: Self) -> float:
        return max(hi - lo for lo, hi in zip(self.min, self.max))

    @property
    # calculate box volume
    def volume(self: Self) -> float:
        sides = [hi - lo for lo, hi in zip(self.min, self.max)]
        return math.prod(sides)

    # return index of dimension with longest side
    def _max_side_idx(self: Self) -> int:
        sides = [hi - lo for lo, hi in zip(self.min, self.max)]
        return sides.index(max(sides))

    # split on longest dimension
    def split_on_longest(self: Self) -> Tuple[Self, Self]:
        idx = self._max_side_idx()
        mid = (self.min[idx] + self.max[idx]) / 2
        # new tuple with updated element at idx
        max1 = self.max[:idx] + (mid,) + self.max[idx + 1 :]
        min2 = self.min[:idx] + (mid,) + self.min[idx + 1 :]

        return BoxN(self.min, Point(max1)), BoxN(Point(min2), self.max)

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
