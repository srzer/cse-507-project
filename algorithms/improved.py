from collections import deque

from typing import Optional, List, Deque
from dreal import And, CheckSatisfiability, Formula, Minimize, Variable

from algorithms import Algorithm
from algorithms.feasible import FeasibleMinBranchAndBound
from box import BoxN
from poly import Rational

# TODO: export these imports in ./box/__init__.py
from box.constraints import build_constraints
from box.feasibility import full_check
from box.split import split_on_longest
from box.type import from_box_model
from poly.bernstein import bernstein_bounds
from poly.type import eval_rational, eval_symbolic


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

        lower_bound = eval_rational(obj, from_box_model(model).center)

        queue: Deque[BoxN] = deque()
        queue.append(init_box)

        while queue:
            box = queue.pop()
            box_constraints = build_constraints(box, vars)

            berstein_min, _ = bernstein_bounds(obj, box)

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
                f_at_mids = eval_rational(obj, from_box_model(improved_model).center)
                if f_at_mids < lower_bound:
                    lower_bound = f_at_mids
                    # print("Updated global lower bound B =", B)

            # 2. If the box is already small enough, skip
            if box.max_side_length <= min_box_size:
                continue

            # FIXME: we probably want the complete feasibility check here,
            # the one with all of the heuristics?
            # instead of just the full check
            if full_check(constr, box_constraints, delta):
                temp_lower_bound = FeasibleMinBranchAndBound(lower_bound)(
                    dim, box, obj, vars, init_constr, min_box_size, delta, err
                )
                if temp_lower_bound and temp_lower_bound < lower_bound:
                    print(
                        "new temp lower bound found with feasible: ", temp_lower_bound
                    )
                    lower_bound = temp_lower_bound
                continue

            b1, b2 = split_on_longest(box)
            queue.append(b1)
            queue.append(b2)

        # print("Final approximate global lower bound B =", B)
        return lower_bound
