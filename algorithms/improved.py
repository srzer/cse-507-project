from collections import deque
from typing import Deque, List

from dreal import And, CheckSatisfiability, Formula, Variable

from box import (
    BoxN,
    BoxSplit,
    build_constraints,
    full_check,
    from_box_model,
)
from objective import (
    Rational,
    ObjectiveBounds,
    eval_rational,
    eval_symbolic,
)

from .either import Either, Left, Right
from .feasible import FeasibleMinBranchAndBound
from .type import Algorithm
from .errors import CONVERGENCE_TOLERANCE, MAX_STAGNANT_ITERS


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
    ) -> Either[str, float]:
        fn_expr = eval_symbolic(obj, vars)
        init_constr = And(build_constraints(init_box, vars), constr)

        model = CheckSatisfiability(init_constr, delta)
        if not model:
            return Left("Model is not satisfiable with given constraints")

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
                print(
                    f"  Converged: no improvement for {MAX_STAGNANT_ITERS} iterations"
                )
                print(f"  Final lower bound: {lower_bound}")
                break

            # log progress periodically
            if iteration_count % 100 == 0:
                stagnant_count = iteration_count - last_improvement_iter
                print(
                    f"  iteration {iteration_count}, queue size: {len(queue)}, current lower bound: {lower_bound}, stagnant: {stagnant_count}"
                )
            box = queue.pop()
            box_constraints = build_constraints(box, vars)

            bounded_min, _ = bounder(obj, box)

            if bounded_min >= lower_bound - err:
                continue
            # 1. Prune by checking if the box can still improve the global lower bound:
            #    Is there any point in this box with constraint satisfied and obj < lower bound - err?

            improved_formula = And(box_constraints, constr, fn_expr < lower_bound - err)
            improved_model = CheckSatisfiability(improved_formula, delta)
            if improved_model is None:
                # no such point ⇒ this box cannot improve lower bound ⇒ discard
                # print("Pruned a box.")
                continue
            else:
                # extract a point from m_improve, and update lower bound
                f_at_mids = eval_rational(obj, from_box_model(improved_model).center)
                # consider bounds "unchanged" if diff < this
                if f_at_mids < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = f_at_mids
                    last_improvement_iter = iteration_count

            # 2. If the box is already small enough, skip
            if box.max_side_length <= min_box_size:
                continue

            # FIXME: we probably want the complete feasibility check here,
            # the one with all of the heuristics?
            # instead of just the full check
            if full_check(constr, box_constraints, delta):
                temp_lower_bound = FeasibleMinBranchAndBound(lower_bound)(
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
                # TODO use a map here instead?
                match temp_lower_bound:
                    case Left(e):
                        return Left(e)
                    case Right(v):
                        if v < lower_bound - CONVERGENCE_TOLERANCE:
                            lower_bound = v
                            last_improvement_iter = iteration_count
                continue

            b1, b2 = splitter(box, obj)
            queue.append(b1)
            queue.append(b2)

        return Right(lower_bound)
