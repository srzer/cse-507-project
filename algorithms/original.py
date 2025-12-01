from collections import deque

from typing import Optional, List
from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from box import BoxN

# TODO: export these imports in ./box/__init__.py
from box.constraints import build_constraints
from box.feasibility import full_check
from box.type import from_box_model
from poly import Rational
from poly.type import eval, eval_symbolic


class GlobalMinBranchAndBound(Algorithm):
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
        print("objective function: ")
        print(fn_expr)

        # feasibility check & initial lower bound
        init_constraints = And(build_constraints(init_box, vars), constr)

        model_box = CheckSatisfiability(init_constraints, delta)
        if model_box is None:
            print("No feasible point in the initial region under the given constraint.")
            return None

        # use the model_box interval midpoints for initial feasible point
        # initial lower bound from feasible point
        lower_bound = eval(obj, from_box_model(model_box).center)
        print("initial feasible point:  ", lower_bound)

        queue = deque()
        queue.append(init_box)
        # branch and bound loop
        while queue:
            box = queue.pop()
            box_constraints = build_constraints(box, vars)

            # TODO: implement bernstein bounds
            # berstein_min = bernstein_bounds_on_box(obj, box)[0]
            # if berstein_min >= lower_bound - err:
            #     print("it works!")
            #     continue
            # 1. Prune by checking if the box can still improve the global bound lower_bound:
            #    Is there any point in this box with constraint satisfied
            #    and f < lower_bound - eps?
            improved_formula = And(box_constraints, constr, fn_expr < lower_bound - err)
            improved_model = CheckSatisfiability(improved_formula, delta)
            if improved_model is None:
                # No such point ⇒ this box cannot improve B ⇒ discard
                # print("Pruned a box.")
                continue
            else:
                # extract a point from m_improve, and update B
                f_at_mids = eval(obj, from_box_model(improved_model).center)

                if f_at_mids < lower_bound:
                    lower_bound = f_at_mids
                    # print("Updated global lower bound B =", B)

            # 2. If the box is already small enough, do a local Minimize on this box
            if box.max_side <= min_box_size:
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
                f_min_approx = eval(obj, mids)

                if f_min_approx < lower_bound:
                    lower_bound = f_min_approx
                #     # print("Updated global lower bound B =", B)
                continue

            # 2bis. If the entire box is within the feasible region, we can also
            #       call Minimize directly on the whole box.

            # fully_feasible = box_is_fully_feasible(box_constraints, constraint, delta_dreal)
            # TODO: implement the whole feasible heuristics function here when ready
            fully_feasible = full_check(constr, box_constraints, delta)
            if fully_feasible:
                # print("Box fully feasible, performing minimize on full box.")
                sol_box = Minimize(fn_expr, box_constraints, delta)
                if not sol_box:
                    return None

                f_min_approx = eval(obj, from_box_model(sol_box).center)

                if f_min_approx < lower_bound:
                    lower_bound = f_min_approx
                    # print("Updated global lower bound B =", B)
                continue

            # 3. Otherwise, the box intersects the feasible region but is not
            #    fully inside it, and it may still contain points with f < B - eps.
            #    We split it and continue the search.
            b1, b2 = box.split_on_longest()
            queue.append(b1)
            queue.append(b2)
            # print("Split box into two children.")

        print("Final approximate global lower bound B =", lower_bound)
        return lower_bound
