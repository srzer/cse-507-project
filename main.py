import time
from collections import deque

from typing import Optional
from dreal import And, CheckSatisfiability, Minimize, Variable

from baseline import baseline_min_dreal
from box import BoxN, Point

# TODO: export these imports in ./box/__init__.py
from box.constraints import build_basic_box, build_constraints
from box.feasibility import full_check
from box.type import max_side, split_on_longest
from poly import Rational
from poly.type import eval, eval_symbolic
from testing.example import ball_constraint_example, rational_objective_example


def global_min_branch_and_bound(
    init_box: BoxN,
    dim: int,
    obj: Rational,
    delta: float = 1e-3,
    min_box_size: float = 0.1,
    err: float = 1e-4,
) -> Optional[float]:
    """
    Global minimization of f(x1,...,xn) under a general constraint using
    a branch-and-bound style algorithm with dReal.

    Parameters
    ----------
    initial_box : BoxND
        Outer search box in n dimensions.
    dim : int
        Dimension (number of variables).
    delta : float
        δ-precision for dReal.
    min_box_size : float
        Stop splitting a box when its longest side is below this threshold.
    eps : float
        Pruning margin: we discard boxes that cannot have f < B - err,
        where B is the current global lower bound.
    """
    assert init_box.dim == dim, "initial_box dimension does not match n."

    # dReal variables x0, x1, ..., x_{n-1}
    vars = [Variable(f"x{i}") for i in range(dim)]

    # symbolic constraint and objective
    constraint = ball_constraint_example(vars)
    fn_expr = eval_symbolic(obj, vars)
    print("objective function: ")
    print(fn_expr)

    # ----- Step 1: Feasibility check + initial lower bound B -----

    init_constraints = And(build_constraints(init_box, vars), constraint)

    model_box = CheckSatisfiability(init_constraints, delta)
    if model_box is None:
        print("No feasible point in the initial region under the given constraint.")
        return None

    # use the model_box interval midpoints for initial feasible point
    intervals = [model_box[x_i] for x_i in vars]
    mids = [iv.mid() for iv in intervals]
    # initial lower bound from feasible point
    lower_bound = eval(obj, Point(tuple(mids)))
    print("initial feasible point:  ", lower_bound)

    # ----- Branch-and-bound main loop -----

    queue = deque()
    queue.append(init_box)

    while queue:
        box = queue.pop()
        box_constraints = build_constraints(box, vars)

        # berstein_min = bernstein_bounds_on_box(obj, box)[0]
        # if berstein_min >= lower_bound - err:
        #     print("it works!")
        #     continue
        # 1. Prune by checking if the box can still improve the global bound B:
        #    Is there any point in this box with constraint satisfied and f < B - eps?
        improve_formula = And(box_constraints, constraint, fn_expr < lower_bound - err)
        m_improve = CheckSatisfiability(improve_formula, delta)
        if m_improve is None:
            # No such point ⇒ this box cannot improve B ⇒ discard
            # print("Pruned a box.")
            continue
        else:
            # extract a point from m_improve, and update B
            intervals = [m_improve[xi] for xi in vars]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_at_mids = eval(obj, Point(tuple(mids)))
            if f_at_mids < lower_bound:
                lower_bound = f_at_mids
                # print("Updated global lower bound B =", B)

        # 2. If the box is already small enough, do a local Minimize on this box
        if max_side(box) <= min_box_size:
            continue
            # NOTE: 3 ways;
            # 1) directly skip
            # 2) minimize with only the box constraint
            # 3) minimize with both box and global constraint
            # NOTE: if we set min_box_size = delta_dreal, then we can directly skip.
            # # print("Box small enough, performing local minimize.")
            # local_constraints = box_constraints
            # # You can choose to also include the global constraint here:
            local_constraints = And(box_constraints, constraint)

            # # dreal
            sol_box = Minimize(fn_expr, local_constraints, delta)
            intervals = [sol_box[xi] for xi in vars]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_min_approx = eval(fn_obj, mids)

            if f_min_approx < lower_bound:
                lower_bound = f_min_approx
            #     # print("Updated global lower bound B =", B)
            continue

        # 2bis. If the entire box is within the feasible region, we can also
        #       call Minimize directly on the whole box.

        # fully_feasible = box_is_fully_feasible(box_constraints, constraint, delta_dreal)
        # TODO: implement the whole feasible heuristics function here when ready
        fully_feasible = full_check(constraint, box_constraints, delta)
        if fully_feasible:
            # print("Box fully feasible, performing minimize on full box.")
            sol_box = Minimize(fn_expr, box_constraints, delta)
            if not sol_box:
                return None

            intervals = [sol_box[xi] for xi in vars]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_min_approx = eval(obj, Point(tuple(mids)))

            if f_min_approx < lower_bound:
                lower_bound = f_min_approx
                # print("Updated global lower bound B =", B)
            continue

        # 3. Otherwise, the box intersects the feasible region but is not
        #    fully inside it, and it may still contain points with f < B - eps.
        #    We split it and continue the search.
        b1, b2 = split_on_longest(box)
        queue.append(b1)
        queue.append(b2)
        # print("Split box into two children.")

    print("Final approximate global lower bound B =", lower_bound)
    return lower_bound


if __name__ == "__main__":
    dim = 3
    init_box = build_basic_box(1, 10, dim)
    fn_obj = rational_objective_example(dim)

    t0 = time.time()
    global_min_branch_and_bound(
        init_box,
        dim,
        fn_obj,
        1e-3,  # δ for dReal
        1e-3,  # stop splitting when max side < min_box_size
        1e-3,  # pruning margin
    )
    t1 = time.time()
    print("Branch-and-bound time:", t1 - t0, "seconds")
    time.sleep(2)

    t2 = time.time()

    baseline_min_dreal(
        init_box,
        dim,
        fn_obj,
        1e-3,  # δ for dReal
    )
    t3 = time.time()
    print("Baseline dReal Minimize time:", t3 - t2, "seconds")
