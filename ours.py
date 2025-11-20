from dataclasses import dataclass
from typing import Tuple, Callable
from collections import deque
from dreal import *
from setting import f_value, f_constraint, get_init_box
import time

# --------- Basic data structures ---------

@dataclass
class Box3D:
    xl: float; xu: float
    yl: float; yu: float
    zl: float; zu: float

def build_box_constraints(x, y, z, box: Box3D):
    """Convert a Box3D into dReal constraints And(...)."""
    return And(box.xl <= x, x <= box.xu,
               box.yl <= y, y <= box.yu,
               box.zl <= z, z <= box.zu)


def max_side(box: Box3D) -> float:
    """Return the maximum side length of the box."""
    return max(box.xu - box.xl, box.yu - box.yl, box.zu - box.zl)


def split_box(box: Box3D) -> Tuple[Box3D, Box3D]:
    """Split the box into two boxes along its longest dimension."""
    dx = box.xu - box.xl
    dy = box.yu - box.yl
    dz = box.zu - box.zl

    if dx >= dy and dx >= dz:
        mid = 0.5 * (box.xl + box.xu)
        b1 = Box3D(box.xl, mid, box.yl, box.yu, box.zl, box.zu)
        b2 = Box3D(mid, box.xu, box.yl, box.yu, box.zl, box.zu)
    elif dy >= dx and dy >= dz:
        mid = 0.5 * (box.yl + box.yu)
        b1 = Box3D(box.xl, box.xu, box.yl, mid, box.zl, box.zu)
        b2 = Box3D(box.xl, box.xu, mid, box.yu, box.zl, box.zu)
    else:
        mid = 0.5 * (box.zl + box.zu)
        b1 = Box3D(box.xl, box.xu, box.yl, box.yu, box.zl, mid)
        b2 = Box3D(box.xl, box.xu, box.yl, box.yu, mid, box.zu)

    return b1, b2


def box_is_fully_feasible(box_constraints,
                          constraint,
                          delta_dreal: float) -> bool:
    """
    Check (using dReal) whether ALL points in this box satisfy `constraint`.

    Idea:
      Box is fully feasible  ⇔  there is NO point in the box that violates constraint
      i.e. UNSAT: box_constraints ∧ ¬constraint
    """
    formula_violate = And(box_constraints, Not(constraint))
    model = CheckSatisfiability(formula_violate, delta_dreal)
    # If UNSAT, there is no violating point → fully feasible (up to δ)
    return model is None


# --------- Main algorithm: branch-and-bound global minimum ---------

def global_min_branch_and_bound(
    initial_box: Box3D,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4
):
    """
    Branch-and-bound style global lower bound for f(x, y, z) under a general constraint.

    Parameters
    ----------
    initial_box : Box3D
        The outer axis-aligned bounding box where we search.

    delta_dreal : float
        δ parameter for dReal (precision of δ-SAT / Minimize).

    min_box_size : float
        If the maximum side length of a box is below this, we stop splitting it.

    eps : float
        Pruning margin. We only keep boxes where it is still δ-possible
        that f < current_B - eps.
    """

    # ---- dReal variables and symbolic expressions ----
    x = Variable("x")
    y = Variable("y")
    z = Variable("z")
    t = Variable("t")  # auxiliary variable for t = f(x, y, z)

    # User-defined global constraint, e.g. sphere, polynomials, etc.
    constraint = f_constraint(x, y, z)

    # Symbolic f(x, y, z)
    f_expr = f_value(x,y,z)

    # ---- Step 1: Feasibility check + initial lower bound B ----

    init_constraints = And(build_box_constraints(x, y, z, initial_box),
                           constraint)

    model = CheckSatisfiability(init_constraints, delta_dreal)
    if model is None:
        print("No feasible point in the initial region under the given constraint.")
        return None

    # Use the midpoint of the model intervals to initialize B
    Ix = model[x]
    Iy = model[y]
    Iz = model[z]

    xm = 0.5 * (Ix.lb() + Ix.ub())
    ym = 0.5 * (Iy.lb() + Iy.ub())
    zm = 0.5 * (Iz.lb() + Iz.ub())

    B = f_value(xm, ym, zm)
    # print("Initial feasible point B =", B)

    # ---- Branch-and-bound main loop ----

    queue = deque()
    queue.append(initial_box)

    while queue:
        box = queue.pop()
        # print(box)
        # Build constraints for this box
        box_constraints = build_box_constraints(x, y, z, box)

        # 1. Pruning by current lower bound:
        #    Check if there exists a point in this box with f < B - eps.
        improve_formula = And(box_constraints, constraint, f_expr < B - eps)
        m_improve = CheckSatisfiability(improve_formula, delta_dreal)
        if m_improve is None:
            # No point in this box can improve the current bound B → prune this box
            # print("cut")
            continue
        # print("could improve")

        # 2. If the whole box satisfies the global constraint
        #    (checked via dReal: box ∧ ¬constraint is UNSAT),
        #    OR the box is already small enough, we run a local Minimize over this box.

        if max_side(box) <= min_box_size:
            # print("block small enough")
            f = f_value(x, y, z)
            new_constraints = box_constraints
            # new_constraints = And(box_constraints, constraint)
            sol_box = Minimize(f, new_constraints, delta_dreal)
            Ix = sol_box[x]
            Iy = sol_box[y]
            Iz = sol_box[z]
            xm = 0.5 * (Ix.lb() + Ix.ub())
            ym = 0.5 * (Iy.lb() + Iy.ub())
            zm = 0.5 * (Iz.lb() + Iz.ub())         
            f_min_approx = f_value(xm, ym, zm)  # approximate lower bound of f in this box

            if f_min_approx < B:
                B = f_min_approx
                # print("Updated global lower bound B =", B)
            continue
        
        fully_feasible = box_is_fully_feasible(
            box_constraints, constraint, delta_dreal
        )
        if fully_feasible:
            # print("fully feasible")
            f = f_value(x, y, z)
            sol_box = Minimize(f, box_constraints, delta_dreal)
            Ix = sol_box[x]
            Iy = sol_box[y]
            Iz = sol_box[z]
            xm = 0.5 * (Ix.lb() + Ix.ub())
            ym = 0.5 * (Iy.lb() + Iy.ub())
            zm = 0.5 * (Iz.lb() + Iz.ub())         
            f_min_approx = f_value(xm, ym, zm)  # approximate lower bound of f in this box

            if f_min_approx < B:
                B = f_min_approx
                # print("Updated global lower bound B =", B)
            continue

        # 3. Otherwise, the box intersects constraint but is not fully inside,
        #    and it might still improve B → split further.
        b1, b2 = split_box(box)
        queue.append(b1)
        queue.append(b2)
        # print("split into", b1, b2)

    print("Final global lower bound (approx) B =", B)
    return B

def baseline_min_dreal(
    initial_box: Box3D,
    delta_dreal: float = 1e-3
):
    """Baseline: directly minimize f over initial_box ∧ f_constraint using dReal."""
    x = Variable("x")
    y = Variable("y")
    z = Variable("z")

    # symbolic f and constraint
    f_expr = f_value(x, y, z)
    constraint = f_constraint(x, y, z)

    # initial region = box ∧ constraint
    region_constraints = And(
        build_box_constraints(x, y, z, initial_box),
        constraint
    )

    # Feasibility check first (optional but safer)
    model = CheckSatisfiability(region_constraints, delta_dreal)
    if model is None:
        print("[baseline] No feasible point in initial_box under f_constraint.")
        return None

    # Direct minimization
    sol_box = Minimize(f_expr, region_constraints, delta_dreal)

    Ix = sol_box[x]
    Iy = sol_box[y]
    Iz = sol_box[z]

    # Use midpoint of the solution box as an approximate minimizer
    xm = 0.5 * (Ix.lb() + Ix.ub())
    ym = 0.5 * (Iy.lb() + Iy.ub())
    zm = 0.5 * (Iz.lb() + Iz.ub())

    f_min_approx = f_value(xm, ym, zm)

    print("[baseline] approximate global minimum f ≈", f_min_approx)

    return f_min_approx


# --------- Example usage ---------

if __name__ == "__main__":
    init_box, min_box_size = get_init_box()
    t0 = time.time()
    global_min_branch_and_bound(
        initial_box=init_box,
        delta_dreal=1e-3,   # δ for dReal
        min_box_size=min_box_size,   # stop splitting when box is this small
        eps=1e-3            # pruning margin
    )
    t1 = time.time()
    print("Branch-and-bound time:", t1 - t0, "seconds")
    baseline_min_dreal(
        initial_box=init_box,
        delta_dreal=1e-3
    )
    t2 = time.time()
    print("Baseline dReal Minimize time:", t2 - t1, "seconds")