from dataclasses import dataclass
from typing import Tuple, List
from collections import deque
from dreal import *
from setting import BoxND, f_value, f_constraint, get_init_box
import time


# --------- Basic helpers for n-D boxes ---------


def build_box_constraints(xs: List[Variable], box: BoxND):
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


# --------- Branch-and-bound global minimization ---------


def global_min_branch_and_bound(
    initial_box: BoxND,
    n: int,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4
):
    """
    Global minimization of f(x1,...,xn) under a general constraint using
    a branch-and-bound style algorithm with dReal.

    Parameters
    ----------
    initial_box : BoxND
        Outer search box in n dimensions.
    n : int
        Dimension (number of variables).
    delta_dreal : float
        δ-precision for dReal.
    min_box_size : float
        Stop splitting a box when its longest side is below this threshold.
    eps : float
        Pruning margin: we discard boxes that cannot have f < B - eps,
        where B is the current global lower bound.
    """
    assert initial_box.dim == n, "initial_box dimension does not match n."

    # dReal variables x0, x1, ..., x_{n-1}
    xs = [Variable(f"x{i}") for i in range(n)]

    # Symbolic constraint and objective
    constraint = f_constraint(*xs)
    f_expr = f_value(*xs)

    # ----- Step 1: Feasibility check + initial lower bound B -----

    init_constraints = And(
        build_box_constraints(xs, initial_box),
        constraint
    )

    model = CheckSatisfiability(init_constraints, delta_dreal)
    if model is None:
        print("No feasible point in the initial region under the given constraint.")
        return None

    # Use the midpoints of model intervals as an initial feasible point
    intervals = [model[xi] for xi in xs]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]

    # Initial lower bound from this feasible point
    B = f_value(*mids)
    # print("Initial feasible point, objective =", B)

    # ----- Branch-and-bound main loop -----

    queue = deque()
    queue.append(initial_box)

    while queue:
        box = queue.pop()

        box_constraints = build_box_constraints(xs, box)

        # 1. Prune by checking if the box can still improve the global bound B:
        #    Is there any point in this box with constraint satisfied and f < B - eps?
        improve_formula = And(box_constraints, constraint, f_expr < B - eps)
        m_improve = CheckSatisfiability(improve_formula, delta_dreal)
        if m_improve is None:
            # No such point ⇒ this box cannot improve B ⇒ discard
            # print("Pruned a box.")
            continue
        else:
            # extract a point from m_improve, and update B
            intervals = [m_improve[xi] for xi in xs]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_at_mids = f_value(*mids)
            if f_at_mids < B:
                B = f_at_mids
                # print("Updated global lower bound B =", B)

        # 2. If the box is already small enough, do a local Minimize on this box
        if max_side(box) <= min_box_size:
            # print("Box small enough, performing local minimize.")
            f_sym = f_value(*xs)
            # local_constraints = box_constraints
            # You can choose to also include the global constraint here:
            local_constraints = And(box_constraints, constraint)

            sol_box = Minimize(f_sym, local_constraints, delta_dreal)

            intervals = [sol_box[xi] for xi in xs]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_min_approx = f_value(*mids)

            if f_min_approx < B:
                B = f_min_approx
                # print("Updated global lower bound B =", B)
            continue

        # 2bis. If the entire box is within the feasible region, we can also
        #       call Minimize directly on the whole box.
        fully_feasible = box_is_fully_feasible(
            box_constraints, constraint, delta_dreal
        )
        if fully_feasible:
            # print("Box fully feasible, performing minimize on full box.")
            f_sym = f_value(*xs)
            sol_box = Minimize(f_sym, box_constraints, delta_dreal)

            intervals = [sol_box[xi] for xi in xs]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_min_approx = f_value(*mids)

            if f_min_approx < B:
                B = f_min_approx
                # print("Updated global lower bound B =", B)
            continue

        # 3. Otherwise, the box intersects the feasible region but is not
        #    fully inside it, and it may still contain points with f < B - eps.
        #    We split it and continue the search.
        b1, b2 = split_box(box)
        queue.append(b1)
        queue.append(b2)
        # print("Split box into two children.")

    print("Final approximate global lower bound B =", B)
    return B


# --------- Baseline: direct dReal Minimize on whole region ---------


def baseline_min_dreal(
    initial_box: BoxND,
    n: int,
    delta_dreal: float = 1e-3
):
    """
    Baseline method:

        Directly call dReal.Minimize on the whole region

            initial_box ∧ f_constraint

        without any branch-and-bound splitting or pruning.
    """
    assert initial_box.dim == n, "initial_box dimension does not match n."

    xs = [Variable(f"x{i}") for i in range(n)]

    f_expr = f_value(*xs)
    constraint = f_constraint(*xs)

    region_constraints = And(
        build_box_constraints(xs, initial_box),
        constraint
    )

    # First check feasibility
    model = CheckSatisfiability(region_constraints, delta_dreal)
    if model is None:
        print("[baseline] No feasible point in initial_box under f_constraint.")
        return None

    # Then perform global Minimize over the region
    sol_box = Minimize(f_expr, region_constraints, delta_dreal)

    intervals = [sol_box[xi] for xi in xs]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]

    f_min_approx = f_value(*mids)

    print("[baseline] approximate global minimum f ≈", f_min_approx)
    return f_min_approx


# --------- Example main ---------


if __name__ == "__main__":
    # Dimension n: you can change this manually or parse from command line
    n = 3

    # get_init_box(n) should be defined in setting.py and return:
    #   BoxND(lows, highs), min_box_size
    init_box, min_box_size = get_init_box(n)

    t0 = time.time()
    global_min_branch_and_bound(
        initial_box=init_box,
        n=n,
        delta_dreal=1e-3,       # δ for dReal
        min_box_size=min_box_size,  # stop splitting when max side < min_box_size
        eps=1e-3                # pruning margin
    )
    t1 = time.time()
    print("Branch-and-bound time:", t1 - t0, "seconds")

    baseline_min_dreal(
        initial_box=init_box,
        n=n,
        delta_dreal=1e-3
    )
    t2 = time.time()
    print("Baseline dReal Minimize time:", t2 - t1, "seconds")
