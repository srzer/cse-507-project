from collections import deque
from typing import Deque, List, Optional

from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from box import BoxN, BoxSplit, FullFeasible, from_box_model
from box.feasibility.check import FullFeasible
from objective import Rational, ObjectiveBounds, eval_rational, eval_symbolic

from .either import Either, Left, Right
from .errors import (
    CONVERGENCE_TOLERANCE,
    MAX_STAGNANT_ITERS,
    ERROR_INFEASIBLE,
    ERROR_NO_GLOBAL_MIN,
)


class GlobalMinBranchAndBound(Algorithm):
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
        Global minimization of f(x1,...,xn) under a general constraint using
        a branch-and-bound style algorithm with dReal.

        Parameters
        ----------
        initial_box : BoxND
            Outer search box in n dimensions.
        dim : int
            Dimension (number of variables).
        delta : float
            δ-precision for dReal.
        min_box_size : float
            Stop splitting a box when its longest side is below this threshold.
        eps : float
            Pruning margin: we discard boxes that cannot have f < B - err,
            where B is the current global lower bound.
        """

        # symbolic constraint and objective
        fn_expr = eval_symbolic(obj, vars)

        # feasibility check & initial lower bound
        init_constraints = And(init_box.build_constraints(vars), constr)

        model_box = CheckSatisfiability(init_constraints, delta)
        if model_box is None:
            return Left(ERROR_INFEASIBLE)

        # use the model_box interval midpoints for initial feasible point
        # initial lower bound from feasible point
        lower_bound = eval_rational(obj, from_box_model(model_box).center)

        queue: Deque[BoxN] = deque()
        queue.append(init_box)

        # branch and bound loop
        iteration_count = 0
        # convergence tracking
        last_improvement_iter = 0

        # give up after this many iterations without improvement

        while queue:
            iteration_count += 1

            # convergence check: give up if no improvement for too long
            if iteration_count - last_improvement_iter > MAX_STAGNANT_ITERS:
                print(
                    f"  converged: no improvement for {MAX_STAGNANT_ITERS} iterations"
                )
                print(f"  final lower bound: {lower_bound}")
                break

            # TODO: implement iteration logging
            # (ideally as monadic state transformer)
            if iteration_count % 100 == 0:
                stagnant_count = iteration_count - last_improvement_iter
                print(
                    f"  iteration {iteration_count}, queue size: {len(queue)}, current lower bound: {lower_bound}, stagnant: {stagnant_count}"
                )

            box: BoxN = queue.pop()
            box_constraints = box.build_constraints(vars)

            # 1. Prune by checking if the box can still improve the global bound lower_bound:
            #    Is there any point in this box with constraint satisfied
            #    and f < lower_bound - eps?
            improved_formula = And(box_constraints, constr, fn_expr < lower_bound - err)
            improved_model = CheckSatisfiability(improved_formula, delta)
            if improved_model is None:
                # no such point ⇒ this box cannot improve lower bound ⇒ discard
                continue
            else:
                # extract a point from improved model, and update lower bound
                f_at_mids = eval_rational(obj, from_box_model(improved_model).center)

                # consider bounds unchanged if diff < convergence tolerance
                if f_at_mids < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = f_at_mids
                    last_improvement_iter = iteration_count

            # 2. If the box is already small enough, do a local Minimize on this box
            if box.max_side_length <= min_box_size:
                continue
                # NOTE: 3 ways;
                # 1) directly skip
                # 2) minimize with only the box constraint
                # 3) minimize with both box and global constraint
                # NOTE: if we set min_box_size = delta_dreal, then we can directly skip.
                # # print("Box small enough, performing local minimize.")
                # local_constraints = box_constraints
                # # You can choose to also include the global constraint here:
                local_constraints = And(box_constraints, constr)

                # # dreal
                sol_box = Minimize(fn_expr, local_constraints, delta)
                intervals = [sol_box[xi] for xi in vars]
                mids = [iv.mid() for iv in intervals]
                f_min_approx = eval_rational(obj, mids)

                if f_min_approx < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = f_min_approx
                    last_improvement_iter = iteration_count
                continue

            # 2bis. If the entire box is within the feasible region, we can also
            #       call Minimize directly on the whole box.

            # TODO: implement the CompleteFeasible heuristics function here when ready
            fully_feasible = FullFeasible()(box, vars, constr, delta)
            if fully_feasible:
                sol_box = Minimize(fn_expr, box_constraints, delta)
                if not sol_box:
                    return Left(ERROR_NO_GLOBAL_MIN)

                f_min_approx = eval_rational(obj, from_box_model(sol_box).center)

                if f_min_approx < lower_bound - CONVERGENCE_TOLERANCE:
                    lower_bound = f_min_approx
                    last_improvement_iter = iteration_count
                continue

            # 3. Otherwise, the box intersects the feasible region but is not
            #    fully inside it, and it may still contain points with obj < lower bound - eps.
            #    We split it and continue the search.
            b1, b2 = splitter(box, obj)
            queue.append(b1)
            queue.append(b2)

        print(f"algorithm completed after {iteration_count} iterations")
        return Right(lower_bound)
