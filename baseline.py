from typing import Optional

from dreal import And, CheckSatisfiability, Minimize, Variable

from box import BoxN, Point
from box.constraints import build_constraints
from poly import Rational
from poly.type import eval, eval_symbolic
from testing.example import ball_constraint_example, rational_objective_example


# baseline algorithm: direct dReal Minimize on whole region
def baseline_min_dreal(
    init_box: BoxN, dim: int, obj: Rational, delta: float = 1e-3
) -> Optional[float]:
    """
    Baseline method:
        Directly call dReal.Minimize on the whole region
            initial_box ∧ f_constraint
        without any branch-and-bound splitting or pruning.
    """
    assert init_box.dim == dim, "initial_box dimension does not match n."

    xs = [Variable(f"x{i}") for i in range(dim)]

    fn_expr = eval_symbolic(obj, xs)
    print(fn_expr)
    constraint = ball_constraint_example(xs)

    region_constraints = And(build_constraints(init_box, xs), constraint)

    # check feasibility
    model = CheckSatisfiability(region_constraints, delta)
    if model is None:
        print("[baseline] No feasible point in initial_box under f_constraint.")
        return None

    # perform global Minimize over the region
    sol_box = Minimize(fn_expr, region_constraints, delta)
    if not sol_box:
        return None

    intervals = [sol_box[xi] for xi in xs]
    mids = [iv.mid() for iv in intervals]
    min_approx = eval(obj, Point(tuple(mids)))

    print("[baseline] approximate global minimum f ≈", min_approx)
    return min_approx
