from dataclasses import dataclass
from typing import Tuple
from dreal import *

@dataclass
class BoxND:
    """Axis-aligned box in n dimensions."""
    lows: list
    highs: list

    def __post_init__(self):
        assert len(self.lows) == len(self.highs), "Dimension mismatch in BoxND."

    @property
    def dim(self):
        return len(self.lows)
# --------- Basic helpers for n-D boxes ---------


def build_box_constraints(xs, box):
    """
    Build dReal constraints for an axis-aligned n-D box:

        ∧_i (low_i <= x_i <= high_i)
    """
    assert len(xs) == box.dim
    cons = []
    for xi, lo, hi in zip(xs, box.lows, box.highs):
        cons.append(lo <= xi)
        cons.append(xi <= hi)
    return And(*cons)


def max_side(box: BoxND) -> float:
    """Return the length of the longest side of the box."""
    return max(hi - lo for lo, hi in zip(box.lows, box.highs))


def split_box(box: BoxND) -> Tuple[BoxND, BoxND]:
    """
    Split a box into two along the longest dimension.

    We choose the dimension k with the largest side length and
    cut it at the midpoint.
    """
    lengths = [hi - lo for lo, hi in zip(box.lows, box.highs)]
    k = max(range(box.dim), key=lambda i: lengths[i])

    mid = 0.5 * (box.lows[k] + box.highs[k])

    lows1 = box.lows.copy()
    highs1 = box.highs.copy()
    lows2 = box.lows.copy()
    highs2 = box.highs.copy()

    highs1[k] = mid
    lows2[k] = mid

    b1 = BoxND(lows1, highs1)
    b2 = BoxND(lows2, highs2)
    return b1, b2


def box_is_fully_feasible(box_constraints,
                          constraint,
                          delta_dreal: float) -> bool:
    """
    Check whether all points in the box satisfy the constraint.

    Idea:
        "The whole box is feasible"  ⇔  there is NO point in the box
        violating the constraint.

    Formally:
        box is fully feasible  ⇔  (box_constraints ∧ ¬constraint) is UNSAT.
    """
    formula_violate = And(box_constraints, Not(constraint))
    model = CheckSatisfiability(formula_violate, delta_dreal)
    # UNSAT ⇒ no violating point ⇒ box fully feasible
    return model is None
