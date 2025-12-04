from dataclasses import dataclass
from typing import Tuple, List, Optional
from collections import deque
from dreal import *
from setting import BoxND, f_constraint, get_init_box, get_poly_terms
from poly_utils import *
from affine_utils import *
from box_utils import *
import time

# --------- Branch-and-bound global minimization ---------

def global_min_branch_and_bound(
    initial_box: BoxND,
    n: int,
    poly_num,
    poly_den,
    vars: Optional[List[Variable]] = None,
    constraint: Optional[Formula] = None,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4,
    split_func = None
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

    if vars is None:
        xs = [Variable(f"x{i}") for i in range(n)]
    else:
        xs = vars

    if constraint is None:
        constraint = f_constraint(*xs)

    # Symbolic constraint and objective
    f_expr = rational_expr_symbolic(poly_num, poly_den, xs)
    # print(f_expr)

    # ----- Step 1: Feasibility check + initial lower bound B -----

    init_constraints = And(
        build_box_constraints(xs, initial_box),
        constraint
    )

    model = CheckSatisfiability(init_constraints, delta_dreal)
    if model is None:
        # print("No feasible point in the initial region under the given constraint.") 
        return None

    # Use the midpoints of model intervals as an initial feasible point
    intervals = [model[xi] for xi in xs]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]

    # Initial lower bound from this feasible point
    B = rational_value_numeric(poly_num, poly_den, mids)
    # print("Initial feasible point, objective =", B)

    # ----- Branch-and-bound main loop -----

    queue = deque()
    queue.append(initial_box)

    while queue:
        box = queue.pop()
        box_constraints = build_box_constraints(xs, box)

        affine_min = affine_bounds_on_box(poly_num, poly_den, box)[0]
        if affine_min >= B-eps:
            continue
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
            f_at_mids = rational_value_numeric(poly_num, poly_den, mids)
            if f_at_mids < B:
                B = f_at_mids
                # print("Updated global lower bound B =", B)

        # 2. If the box is already small enough, do a local Minimize on this box
        if max_side(box) <= min_box_size:
            continue

        # 2bis. If the entire box is within the feasible region, we can also
        #       call Minimize directly on the whole box.
        fully_feasible = box_is_fully_feasible(
            box_constraints, constraint, delta_dreal
        )
        if fully_feasible:
            # print("Box fully feasible, performing minimize on full box.")
            sol_box = Minimize(f_expr, box_constraints, delta_dreal)
            if sol_box:
                intervals = [sol_box[xi] for xi in xs]
                mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
                f_min_approx = rational_value_numeric(poly_num, poly_den, mids)
                if f_min_approx < B:
                    B = f_min_approx
            continue

        # 3. Otherwise, the box intersects the feasible region but is not
        #    fully inside it, and it may still contain points with f < B - eps.
        #    We split it and continue the search.
        b1, b2 = split_func(poly_num, poly_den, box)
        queue.append(b1)
        queue.append(b2)
        # print("Split box into two children.")

    # print("Final approximate global lower bound B =", B) 
    return B


# --------- Improved ----------

def feasible_min_branch_and_bound(
    initial_box: BoxND,
    n: int,
    poly_num,
    poly_den,
    vars: Optional[List[Variable]] = None,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4,
    initial_B: float = 1e9,
    split_func = None
):
    assert initial_box.dim == n, "initial_box dimension does not match n."

    # dReal variables x0, x1, ..., x_{n-1}
    if vars is None:
        xs = [Variable(f"x{i}") for i in range(n)]
    else:
        xs = vars

    f_expr = rational_expr_symbolic(poly_num, poly_den, xs)

    B = initial_B

    # ----- Branch-and-bound main loop -----

    queue = deque()
    queue.append(initial_box)

    while queue:
        box = queue.pop()
        
        affine_min = affine_bounds_on_box(poly_num, poly_den, box)[0]
        if affine_min >= B-eps:
            continue
        corners = []
        def generate_corners(idx, current):
            if idx == box.dim:
                corners.append(current.copy())
                return
            current.append(box.lows[idx])
            generate_corners(idx + 1, current)
            current.pop()
            current.append(box.highs[idx])
            generate_corners(idx + 1, current)
            current.pop()
        generate_corners(0, [])
        for corner in corners:
            f_at_corner = rational_value_numeric(poly_num, poly_den, corner)
            if f_at_corner < B:
                B = f_at_corner
        if max_side(box) <= min_box_size:
            intervals = [Interval(lo, hi) for lo, hi in zip(box.lows, box.highs)]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            f_at_mids = rational_value_numeric(poly_num, poly_den, mids)
            if f_at_mids < B:
                B = f_at_mids
            continue

        b1, b2 = split_func(poly_num, poly_den, box)
        queue.append(b1)
        queue.append(b2)
    return B

def improved_global_min_branch_and_bound(
    initial_box: BoxND,
    n: int,
    poly_num,
    poly_den,
    vars: Optional[List[Variable]] = None,
    constraint: Optional[Formula] = None,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4,
    split_func = None,
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

    if vars is None:
        xs = [Variable(f"x{i}") for i in range(n)]
    else:
        xs = vars

    if constraint is None:
        constraint = f_constraint(*xs)

    f_expr = rational_expr_symbolic(poly_num, poly_den, xs)
    # print(f_expr) 

    init_constraints = And(
        build_box_constraints(xs, initial_box),
        constraint
    )

    model = CheckSatisfiability(init_constraints, delta_dreal)
    if model is None:
        # print("No feasible point in the initial region under the given constraint.") 
        return None

    # Use the midpoints of model intervals as an initial feasible point
    intervals = [model[xi] for xi in xs]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]

    # Initial lower bound from this feasible point
    B = rational_value_numeric(poly_num, poly_den, mids)
    # print("Initial feasible point, objective =", B)

    # ----- Branch-and-bound main loop -----

    queue = deque()
    queue.append(initial_box)

    while queue:
        box = queue.pop()
        box_constraints = build_box_constraints(xs, box)

        affine_min = affine_bounds_on_box(poly_num, poly_den, box)[0]
        if affine_min >= B-eps:
            continue
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
            f_at_mids = rational_value_numeric(poly_num, poly_den, mids)
            if f_at_mids < B:
                B = f_at_mids
                # print("Updated global lower bound B =", B)

        # 2. If the box is already small enough, skip
        if max_side(box) <= min_box_size:
            continue
        
        feasible_flag = box_is_fully_feasible(
            box_constraints, constraint, delta_dreal
        )
        if feasible_flag:
            B_temp = feasible_min_branch_and_bound(
                box,
                n,
                poly_num,
                poly_den,
                vars=xs,
                delta_dreal=delta_dreal,
                min_box_size=min_box_size,
                eps=eps,
                initial_B=B,
                split_func=split_func
            )
            if B_temp < B:
                B = B_temp
            continue
            
        b1, b2 = split_func(poly_num, poly_den, box)
        queue.append(b1)
        queue.append(b2)

    # print("Final approximate global lower bound B =", B) 
    return B

# --------- Baseline: direct dReal Minimize on whole region ---------

def baseline_min_dreal(
    initial_box: BoxND,
    n: int,
    poly_num,
    poly_den,
    vars: Optional[List[Variable]] = None,
    constraint: Optional[Formula] = None,
    delta_dreal: float = 1e-3
):
    """
    Baseline method:

        Directly call dReal.Minimize on the whole region

            initial_box ∧ f_constraint

        without any branch-and-bound splitting or pruning.
    """
    assert initial_box.dim == n, "initial_box dimension does not match n."

    if vars is None:
        xs = [Variable(f"x{i}") for i in range(n)]
    else:
        xs = vars

    if constraint is None:
        constraint = f_constraint(*xs)

    f_expr = rational_expr_symbolic(poly_num, poly_den, xs)
    # print(f_expr)

    region_constraints = And(
        build_box_constraints(xs, initial_box),
        constraint
    )

    model = CheckSatisfiability(region_constraints, delta_dreal)
    if model is None:
        # print("[baseline] No feasible point in initial_box under f_constraint.")
        return None

    sol_box = Minimize(f_expr, region_constraints, delta_dreal)
    if not sol_box:
        return None

    intervals = [sol_box[xi] for xi in xs]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
    f_min_approx = rational_value_numeric(poly_num, poly_den, mids)

    # print("[baseline] approximate global minimum f ≈", f_min_approx)
    return f_min_approx


# --------- Example main ---------

if __name__ == "__main__":
    # Dimension n: you can change this manually or parse from command line
    n = 2

    # get_init_box(n) should be defined in setting.py and return:
    #   BoxND(lows, highs), min_box_size
    init_box = get_init_box(n)
    poly_num_terms, poly_den_terms = get_poly_terms(n)
    poly_num = poly_from_terms(poly_num_terms)
    poly_den = poly_from_terms(poly_den_terms)

    t0 = time.time()
    global_min_branch_and_bound(
        initial_box=init_box,
        n=n,
        poly_num=poly_num,
        poly_den=poly_den,
        delta_dreal=1e-3,       # δ for dReal
        min_box_size=1e-3,  # stop splitting when max side < min_box_size
        eps=1e-3,                # pruning margin
        split_func=split_box_length
    )
    t1 = time.time()
    print("Branch-and-bound time:", t1 - t0, "seconds")
    
    time.sleep(2)

    t0 = time.time()
    global_min_branch_and_bound(
        initial_box=init_box,
        n=n,
        poly_num=poly_num,
        poly_den=poly_den,
        delta_dreal=1e-3,       # δ for dReal
        min_box_size=1e-3,  # stop splitting when max side < min_box_size
        eps=1e-3,                # pruning margin
        split_func=split_box_gradient
    )
    t1 = time.time()
    print("Branch-and-bound (heuristic) time:", t1 - t0, "seconds")
    
    time.sleep(2)

    t0 = time.time()
    improved_global_min_branch_and_bound(
        initial_box=init_box,
        n=n,
        poly_num=poly_num,
        poly_den=poly_den,
        delta_dreal=1e-3,       # δ for dReal
        min_box_size=1e-3,  # stop splitting when max side < min_box_size
        eps=1e-3,                # pruning margin
        split_func=split_box_length
    )
    t1 = time.time()
    print("Improved Branch-and-bound time:", t1 - t0, "seconds")
    
    time.sleep(2)

    t0 = time.time()
    improved_global_min_branch_and_bound(
        initial_box=init_box,
        n=n,
        poly_num=poly_num,
        poly_den=poly_den,
        delta_dreal=1e-3,       # δ for dReal
        min_box_size=1e-3,  # stop splitting when max side < min_box_size
        eps=1e-3,                # pruning margin
        split_func=split_box_gradient
    )
    t1 = time.time()
    print("Improved Branch-and-bound (heuristic) time:", t1 - t0, "seconds")

    t4 = time.time()
    baseline_min_dreal(
        initial_box=init_box,
        n=n,
        poly_num=poly_num,
        poly_den=poly_den,
        delta_dreal=1e-3
    )
    t5 = time.time()
    print("Baseline dReal Minimize time:", t5 - t4, "seconds")