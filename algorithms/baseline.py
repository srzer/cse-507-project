from typing import List

from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm

# FIXME: export the dot imports in the respective __init__ files
from algorithms.either import Either, Left, Right
from box import BoxN
from box.constraints import build_constraints
from box.type import from_box_model
from poly import Rational
from poly.type import eval_rational, eval_symbolic
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
    ) -> Either[str, float]:
        """
        Directly call dReal.Minimize on the whole region
            initial_box ∧ f_constraint
        without any branch-and-bound splitting or pruning.
        """
        fn_expr = eval_symbolic(obj, vars)
        constraint = ball_constraint_example(vars)

        region_constraints = And(build_constraints(init_box, vars), constraint)

        # check feasibility
        model = CheckSatisfiability(region_constraints, delta)
        if model is None:
            return Left(
                "No feasible point in the inital box under the given constraint"
            )

        # perform global Minimize over the region
        sol_box = Minimize(fn_expr, region_constraints, delta)
        if not sol_box:
            return Left("No global min found over given region")

        min_approx = eval_rational(obj, from_box_model(sol_box).center)

        # print("[baseline] approximate global minimum f ≈", min_approx)
        return Right(min_approx)
