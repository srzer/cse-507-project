from collections import deque

from typing import Optional, List
from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from algorithms.feasible import FeasibleMinBranchAndBound
from box import BoxN, Point

# TODO: export these imports in ./box/__init__.py
from box.constraints import build_constraints
from box.feasibility import full_check
from box.type import from_box_model
from poly import Rational
from poly.type import eval, eval_symbolic


class ImprovedGlobalMinBranchAndBound(Algorithm):
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
        fn_expr = eval_symbolic(obj, vars)
        init_constr = And(build_constraints(init_box, vars), constr)

        model = CheckSatisfiability(init_constr, delta)
        if not model:
            return None

        lower_bound = eval(obj, from_box_model(model).center)

        queue = deque()
        queue.append(init_box)

        while queue:
            box = queue.pop()
            box_constraints = build_constraints(box, vars)

            # TODO: implement bernstein!
            berstein_min = bernstein_bounds_on_box(poly_num, poly_den, box)[0]

            if berstein_min >= lower_bound - err:
                continue
            # 1. Prune by checking if the box can still improve the global bound B:
            #    Is there any point in this box with constraint satisfied and f < B - eps?

            improved_formula = And(box_constraints, constr, fn_expr < lower_bound - err)
            improved_model = CheckSatisfiability(improved_formula, delta)
            if improved_model is None:
                # No such point ⇒ this box cannot improve B ⇒ discard
                # print("Pruned a box.")
                continue
            else:
                # extract a point from m_improve, and update lower bound
                f_at_mids = eval(obj, from_box_model(improved_model).center)
                if f_at_mids < lower_bound:
                    lower_bound = f_at_mids
                    # print("Updated global lower bound B =", B)

            # 2. If the box is already small enough, skip
            if box.max_side <= min_box_size:
                continue

            # FIXME: we probably want the complete feasibility check here,
            # the one with all of the heuristics?
            # instead of just the full check
            if full_check(constr, box_constraints, delta):
                feasible_min_bb = FeasibleMinBranchAndBound()

                feasible_min_bb.initial_lower_bound = lower_bound

                temp_lower_bound = feasible_min_bb(
                    dim, box, obj, vars, init_constr, min_box_size, delta, err
                )
                if temp_lower_bound < lower_bound:
                    lower_bound = temp_lower_bound
                continue

            b1, b2 = box.split_on_longest()
            queue.append(b1)
            queue.append(b2)

        # print("Final approximate global lower bound B =", B)
        return lower_bound
