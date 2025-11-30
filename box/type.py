import math
from dataclasses import dataclass
from typing import NewType, Tuple

Point = NewType("Point", Tuple[float, ...])


# immutable
@dataclass(frozen=True)
class BoxN:
    min: Point
    max: Point

    def __post_init__(self):
        assert len(self.min) == len(self.max), "Box dimension mismatch."

    @property
    def dim(self):
        return len(self.min)


# def dim(box: Box) -> int:
#     assert len(box.min) == len(box.max)
#     return len(box.min)


# return length of longest side
def max_side(box: BoxN) -> float:
    return max(hi - lo for lo, hi in zip(box.min, box.max))


# return index of dimension with longest side
def max_side_idx(box: BoxN) -> int:
    sides = [hi - lo for lo, hi in zip(box.min, box.max)]
    return sides.index(max(sides))


# split on longest dimension
def split_on_longest(box: BoxN) -> Tuple[BoxN, BoxN]:
    idx = max_side_idx(box)
    mid = (box.min[idx] + box.max[idx]) / 2
    # new tuple with updated element at idx
    max1 = box.max[:idx] + (mid,) + box.max[idx + 1 :]
    min2 = box.min[:idx] + (mid,) + box.min[idx + 1 :]

    return BoxN(box.min, Point(max1)), BoxN(Point(min2), box.max)


# calculate box volume
def volume(box: BoxN) -> float:
    sides = [hi - lo for lo, hi in zip(box.min, box.max)]
    return math.prod(sides)


# get center point of box
def center(box: BoxN) -> Point:
    return Point(tuple((hi + lo) / 2 for lo, hi in zip(box.min, box.max)))


# check if box containts a point
def contains(box: BoxN, point: Point) -> bool:
    if len(point) != box.dim:
        raise ValueError(
            f"Point dimension {len(point)} does not match Box dimension {box.dim}"
        )
    return all(lo <= pt <= hi for pt, lo, hi in zip(point, box.min, box.max))
