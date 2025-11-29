from collections import deque
from dreal import *
from box import Box3, build_constraints, max_side, split

def box_is_fully_feasible_test(box_constraints, constraint, delta_dreal: float) -> bool:
    formula_violate = And(box_constraints, Not(constraint))
    model = CheckSatisfiability(formula_violate, delta_dreal)
    return model is None

# ============================================================
#  Main Solver (Exact Match to ours.py logic)
# ============================================================

def global_min_branch_and_bound(
    objective_func,
    constraint_func,
    initial_box: Box3D,
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4,
):
    x, y, z = Variable("x"), Variable("y"), Variable("z")
    
    # 1. Setup Expressions
    constraint = constraint_func(x, y, z)
    f_expr = objective_func(x, y, z)

    # 2. Initial Feasibility Check
    init_constraints = And(build_constraints(initial_box, [x, y, z]), constraint)
    model = CheckSatisfiability(init_constraints, delta_dreal)
    if model is None:
        return None

    # 3. Initialize B (The "Incumbent")
    # We execute this exactly like ours.py: Calculate center, evaluate, set B.
    Ix, Iy, Iz = model[x], model[y], model[z]
    xm = 0.5 * (Ix.lb() + Ix.ub())
    ym = 0.5 * (Iy.lb() + Iy.ub())
    zm = 0.5 * (Iz.lb() + Iz.ub())
    
    # Direct float conversion (No try/except masking)
    B = float(objective_func(xm, ym, zm))

    # 4. Main Loop
    queue = deque([initial_box])

    while queue:
        box = queue.pop()
        box_constraints = build_constraints(box, [x, y, z])

        # Pruning: Is best case in box < current B?
        improve_formula = And(box_constraints, constraint, f_expr < B - eps)
        if CheckSatisfiability(improve_formula, delta_dreal) is None:
            continue

        # Check if we should run local minimization
        run_local = False
        if max_side(box) <= min_box_size:
            run_local = True
        elif box_is_fully_feasible_test(box_constraints, constraint, delta_dreal):
            run_local = True

        if run_local:
            # Local Minimize
            sol_box = Minimize(f_expr, And(box_constraints, constraint), delta_dreal)
            if sol_box:
                Ix, Iy, Iz = sol_box[x], sol_box[y], sol_box[z]
                xm = 0.5 * (Ix.lb() + Ix.ub())
                ym = 0.5 * (Iy.lb() + Iy.ub())
                zm = 0.5 * (Iz.lb() + Iz.ub())
                
                # Direct update check
                val = float(objective_func(xm, ym, zm))
                if val < B:
                    B = val
            continue

        # Split and continue
        b1, b2 = split(box)
        queue.append(b1)
        queue.append(b2)

    return B

def baseline_min_dreal(objective_func, constraint_func, initial_box: Box3D, delta_dreal: float = 1e-3):
    x, y, z = Variable("x"), Variable("y"), Variable("z")
    f_expr = objective_func(x, y, z)
    constraint = constraint_func(x, y, z)
    region = And(build_constraints(initial_box, [x, y, z]), constraint)
    
    if CheckSatisfiability(region, delta_dreal) is None:
        return None
    sol_box = Minimize(f_expr, region, delta_dreal)
    if not sol_box: return None
    
    Ix, Iy, Iz = sol_box[x], sol_box[y], sol_box[z]
    xm = 0.5 * (Ix.lb() + Ix.ub())
    ym = 0.5 * (Iy.lb() + Iy.ub())
    zm = 0.5 * (Iz.lb() + Iz.ub())
    return float(objective_func(xm, ym, zm))
