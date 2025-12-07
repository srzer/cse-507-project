from collections import deque
from typing import Deque, List, Tuple

from dreal import And, CheckSatisfiability, Formula, Variable
from returns.result import Result, Failure, Success

from box import (
    BoxN,
    BoxSplit,
    FullFeasible,
    from_box_model,
)
from objective import (
    ObjectiveBounds,
    Rational,
    eval_rational,
    eval_symbolic,
)

from .errors import CONVERGENCE_TOLERANCE, MAX_STAGNANT_ITERS, ERROR_INFEASIBLE
from .feasible import FeasibleMinBranchAndBound
from .log import LogEntry
from .type import Algorithm


class ImprovedGlobalMinBranchAndBound(Algorithm):
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
        logs: List[LogEntry] = []
        fn_expr = eval_symbolic(obj, vars)
        init_constr = And(init_box.build_constraints(vars), constr)

        model = CheckSatisfiability(init_constr, delta)
        if not model:
            return Failure(ERROR_INFEASIBLE)

        lower_bound = eval_rational(obj, from_box_model(model).center)

        queue: Deque[BoxN] = deque()
        queue.append(init_box)

        # convergence tracking
        iteration_count = 0
        last_improvement_iter = 0

        while queue:
            iteration_count += 1

            # convergence check: give up if no improvement for too long
            if iteration_count - last_improvement_iter > MAX_STAGNANT_ITERS:
                logs.append(
                    {
                        "iteration": iteration_count,
                        "action": "converged",
                        "box": "",
                        "bound": lower_bound,
                        "volume": 0.0,
                        "notes": f"no improvement for {MAX_STAGNANT_ITERS} iters",
                    }
                )
                break

            box = queue.pop()
            box_constraints = box.build_constraints(vars)

            bounded_min, _ = bounder(obj, box)

            if bounded_min >= lower_bound - err:
                logs.append(
                    {
                        "iteration": iteration_count,
                        "action": "prune",
                        "box": str(box),
                        "bound": lower_bound,
                        "volume": box.volume,
                        "notes": "box cannot improve lower bound (bounded_min)",
                    }
                )
                continue
            # 1. Prune by checking if the box can still improve the global lower bound:
            #    Is there any point in this box with constraint satisfied and obj < lower bound - err?

            improved_formula = And(box_constraints, constr, fn_expr < lower_bound - err)
            improved_model = CheckSatisfiability(improved_formula, delta)
            if improved_model is None:
                # no such point ⇒ this box cannot improve lower bound ⇒ discard
                logs.append(
                    {
                        "iteration": iteration_count,
                        "action": "prune",
                        "box": str(box),
                        "bound": lower_bound,
                        "volume": box.volume,
                        "notes": "box cannot improve lower bound (improved_model)",
                    }
                )
                continue
            else:
                # extract a point from m_improve, and update lower bound
                f_at_mids = eval_rational(obj, from_box_model(improved_model).center)
                # consider bounds "unchanged" if diff < this
                if f_at_mids < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = f_at_mids
                    last_improvement_iter = iteration_count
                    logs.append(
                        {
                            "iteration": iteration_count,
                            "action": "update_bound",
                            "box": str(box),
                            "bound": lower_bound,
                            "volume": box.volume,
                            "notes": f"new lower bound: {lower_bound}",
                        }
                    )

            # 2. If the box is already small enough, skip
            if box.max_side_length <= min_box_size:
                logs.append(
                    {
                        "iteration": iteration_count,
                        "action": "prune",
                        "box": str(box),
                        "bound": lower_bound,
                        "volume": box.volume,
                        "notes": f"box is small enough (max side: {box.max_side_length})",
                    }
                )
                continue

            # FIXME: we probably want the complete feasibility check here,
            # the one with all of the heuristics?
            # instead of just the full check
            if FullFeasible()(box, vars, constr, delta):
                result = FeasibleMinBranchAndBound(lower_bound)(
                    dim,
                    box,
                    obj,
                    vars,
                    init_constr,
                    splitter,
                    bounder,
                    min_box_size,
                    delta,
                    err,
                )
                if isinstance(result, Failure):
                    return result

                v, sub_logs = result.unwrap()
                logs.extend(sub_logs)
                if v < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = v
                    last_improvement_iter = iteration_count
                    logs.append(
                        {
                            "iteration": iteration_count,
                            "action": "minimize_feasible",
                            "box": str(box),
                            "bound": lower_bound,
                            "volume": box.volume,
                            "notes": f"box is fully feasible, new lower bound: {lower_bound}",
                        }
                    )

                continue

            b1, b2 = splitter(box, obj)
            queue.append(b1)
            queue.append(b2)
            logs.append(
                {
                    "iteration": iteration_count,
                    "action": "split",
                    "box": str(box),
                    "bound": lower_bound,
                    "volume": box.volume,
                    "notes": f"split into {b1} and {b2}",
                }
            )

        return Success((lower_bound, logs))
