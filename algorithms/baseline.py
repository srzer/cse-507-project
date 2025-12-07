from typing import List, Tuple

from dreal import And, CheckSatisfiability, Formula, Minimize, Variable
from returns.result import Result, Failure, Success

from algorithms import Algorithm
from box import BoxN, BoxSplit, from_box_model
from objective import Rational, ObjectiveBounds, eval_rational, eval_symbolic

from .errors import ERROR_INFEASIBLE, ERROR_NO_GLOBAL_MIN
from .log import LogEntry


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
    ) -> Result[Tuple[float, List[LogEntry]], str]:
        """
        Directly call dReal.Minimize on the whole region
            initial_box ∧ f_constraint
        without any branch-and-bound splitting or pruning.
        """
        fn_expr = eval_symbolic(obj, vars)

        region_constraints = And(init_box.build_constraints(vars), constr)

        # check feasibility
        model = CheckSatisfiability(region_constraints, delta)
        if model is None:
            return Failure(ERROR_INFEASIBLE)

        # perform global dReal minimize over the region
        sol_box = Minimize(fn_expr, region_constraints, delta)
        if not sol_box:
            return Failure(ERROR_NO_GLOBAL_MIN)

        # return function evaluated at solution box center
        return Success((eval_rational(obj, from_box_model(sol_box).center), []))
