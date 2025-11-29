"""Functional box types for n-dimensional and 3D optimization."""

from typing import List, Tuple, NamedTuple
from dreal import *
import numpy as np


class Box(NamedTuple):
    lower: Tuple[float, ...]
    upper: Tuple[float, ...]


class Box3(Box):
    lower: Tuple[float, float, float]  # (xl, yl, zl)
    upper: Tuple[float, float, float]  # (xu, yu, zu)


# get box dimnesions
def dim(box: Box) -> int:
    assert len(box.lower) == len(box.upper)
    return len(box.lower)


# cast Box3 to Box
def to_box(box3: Box3) -> Box:
    return Box(box3.lower, box3.upper)


# attempt coerce Box to Box3
def to_box3(box: Box) -> Box3:
    if len(box.lower) != 3:
        raise ValueError(f"Cannot convert {len(box.lower)}-dimensional box to Box3")
    else:
        # assert len(box.lower) == len(box.upper) == 3
        return Box3(box.lower, box.upper)


# return length of longest side
def max_side(box: Box) -> float:
    return max(u - l for l, u in zip(box.lower, box.upper))


# return index of dimension with longest side
def max_side_idx(box: Box) -> int:
    sides = [u - l for l, u in zip(box.lower, box.upper)]
    return sides.index(max(sides))


# split box on longest dimension
def split(box: Box) -> Tuple[Box, Box]:
    idx = max_side_idx(box)
    mid = (box.lower[idx] + box.upper[idx]) / 2

    lower1 = box.lower
    upper1 = box.upper[:idx] + (mid,) + box.upper[idx + 1 :]
    lower2 = box.lower[:idx] + (mid,) + box.lower[idx + 1 :]
    upper2 = box.upper

    # return the same type as input
    if isinstance(box, Box3):
        return Box3(lower1, upper1), Box3(lower2, upper2)
    else:
        return Box(lower1, upper1), Box(lower2, upper2)


# calculate box volume
def volume(box: Box) -> float:
    sides = [u - l for l, u in zip(box.lower, box.upper)]
    return np.prod(sides)


# get center point of box
def center(box: Box) -> Tuple[float, ...]:
    return tuple((u + l) / 2 for l, u in zip(box.lower, box.upper))


# check if box containts a point
def contains(box: Box, point: Tuple[float, ...]) -> bool:
    lows, highs = box.lower, box.upper

    if len(point) != len(lows):
        raise ValueError(
            f"Point dimension {len(point)} doesn't match box dimension {len(lows)}"
        )

    return all(l <= p <= h for p, l, h in zip(point, lows, highs))


# build dReal constraints for box
def constraints(box: Box, variables: List[Variable]) -> Formula:
    lows, highs = box.lower, box.upper

    if len(variables) != len(lows):
        raise ValueError(f"Expected {len(lows)} variables, got {len(variables)}")

    constraints = []
    for var, lo, hi in zip(variables, lows, highs):
        constraints.append(lo <= var)
        constraints.append(var <= hi)

    return And(*constraints) if constraints else Formula(True)


# test if entire box satisfies given constraint formula phi using a grid of test points
def fully_feasible(
    box: Box, constraint: Formula, variables: List[Variable], epsilon: float = 0.1
) -> bool:
    lows, highs = box.lower, box.upper
    box_dim = len(lows)
    # FIXME: simple choice, where did this come from?
    samples_per_dim = max(2, int(10 / box_dim))
    # Generate grid points
    grid_points = [
        np.linspace(lows[i], highs[i], samples_per_dim) for i in range(box_dim)
    ]

    # Test all grid points
    for point_idx in np.ndindex(*[samples_per_dim] * box_dim):
        test_point = [grid_points[i][point_idx[i]] for i in range(box_dim)]
        # Create substitution map
        subst = {var: val for var, val in zip(variables, test_point)}
        # Check if point satisfies constraint
        if not constraint.Substitute(subst).Evaluate():
            return False

    return True
