from typing import List

from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from box import BoxN, BoxSplit, build_constraints, from_box_model
from objective import Rational, ObjectiveBounds, eval_rational, eval_symbolic

from .either import Either, Left, Right
from .errors import ERROR_INFEASIBLE, ERROR_NO_GLOBAL_MIN


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
        splitter: BoxSplit,
        bounder: ObjectiveBounds,
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

        region_constraints = And(build_constraints(init_box, vars), constr)

        # check feasibility
        model = CheckSatisfiability(region_constraints, delta)
        if model is None:
            return Left(ERROR_INFEASIBLE)

        # perform global dReal minimize over the region
        sol_box = Minimize(fn_expr, region_constraints, delta)
        if not sol_box:
            return Left(ERROR_NO_GLOBAL_MIN)

        # return function evaluated at solution box center
        return Right(eval_rational(obj, from_box_model(sol_box).center))
