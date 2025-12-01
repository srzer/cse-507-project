from typing import List, Optional

from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from box import BoxN
from poly import Rational

# FIXME: export the dot imports in the respective __init__ files
from box.constraints import build_constraints
from box.type import from_box_model
from poly.type import eval, eval_symbolic
from testing.example import ball_constraint_example


# NOTE: doesn't take in some of the params given by iface
# baseline algorithm: direct dReal Minimize on whole region
class BaselineMin(Algorithm):
    def _run(
        self,
        dim: int,
        init_box: BoxN,
        obj: Rational,
        vars: List[Variable],
        constr: Formula,
        min_box_size: float,
        delta: float,
        err: float,
    ) -> Optional[float]:
        """
        Directly call dReal.Minimize on the whole region
            initial_box ∧ f_constraint
        without any branch-and-bound splitting or pruning.
        """
        fn_expr = eval_symbolic(obj, vars)
        print(fn_expr)
        constraint = ball_constraint_example(vars)

        region_constraints = And(build_constraints(init_box, vars), constraint)

        # check feasibility
        model = CheckSatisfiability(region_constraints, delta)
        if model is None:
            print("[baseline] No feasible point in initial_box under f_constraint.")
            return None

        # perform global Minimize over the region
        sol_box = Minimize(fn_expr, region_constraints, delta)
        if not sol_box:
            return None

        min_approx = eval(obj, from_box_model(sol_box).center)

        print("[baseline] approximate global minimum f ≈", min_approx)
        return min_approx
