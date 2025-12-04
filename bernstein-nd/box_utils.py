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


def split_box_length(poly_num, poly_den, box: BoxND) -> Tuple[BoxND, BoxND]:
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


def split_box_gradient(poly_num, poly_den, box: "BoxND") -> Tuple["BoxND", "BoxND"]:
    def eval_rational_grad(poly_num, poly_den, point):
        n = len(point)
        num_val = 0.0
        den_val = 0.0
        grad_num = [0.0] * n
        grad_den = [0.0] * n

        for mon, coeff in poly_num.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            num_val += term

            for k in range(n):
                e_k = mon[k]
                if e_k == 0:
                    continue
                term_grad = coeff * e_k * (point[k] ** (e_k - 1))
                for j, p_j in enumerate(point):
                    if j == k:
                        continue
                    term_grad *= p_j ** mon[j]
                grad_num[k] += term_grad

        for mon, coeff in poly_den.items():
            term = coeff
            for j, p_j in enumerate(point):
                term *= p_j ** mon[j]
            den_val += term

            for k in range(n):
                e_k = mon[k]
                if e_k == 0:
                    continue
                term_grad = coeff * e_k * (point[k] ** (e_k - 1))
                for j, p_j in enumerate(point):
                    if j == k:
                        continue
                    term_grad *= p_j ** mon[j]
                grad_den[k] += term_grad
        den_val = max(den_val, 1e-12)  # prevent division by zero
        grad_f = []
        for k in range(n):
            grad_f_k = (grad_num[k] * den_val - num_val * grad_den[k]) / (den_val ** 2)
            grad_f.append(grad_f_k)

        return grad_f

    mid_point = [(lo + hi) * 0.5 for lo, hi in zip(box.lows, box.highs)]

    grad = eval_rational_grad(poly_num, poly_den, mid_point)

    if any(abs(g) > 0.0 for g in grad):
        k = max(range(box.dim), key=lambda i: abs(grad[i]))
    else:
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
