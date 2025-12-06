from collections import deque

from typing import Optional, List, Deque
from dreal import And, CheckSatisfiability, Formula, Minimize, Variable, Interval

from algorithms import Algorithm
from box import BoxN, Point

# TODO: export these imports in ./box/__init__.py
from box.constraints import build_constraints
from box.feasibility import full_check
from box.split import split_on_longest
from poly import Rational
from poly.type import eval_rational, eval_symbolic
from poly.bernstein import bernstein_bounds

# FIXME: why do we have separate algorithm implementations?
# shouldn't we just have all the heurstic pieces (or other components)
# just be passables to plug and play?
# ... i suppose we can do this for now


class FeasibleMinBranchAndBound(Algorithm):
    def __init__(self, initial_lower_bound: float = float('inf')):
        self.initial_lower_bound = initial_lower_bound

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
        # fn_expr = eval_symbolic(obj, vars)
        lower_bound = self.initial_lower_bound

        queue: Deque[BoxN] = deque()
        queue.append(init_box)

        while queue:
            box: BoxN = queue.pop()

            berstein_min, _ = bernstein_bounds(obj, box)

            if berstein_min >= lower_bound - err:
                continue
            # # 1. Prune by checking if the box can still improve the global bound B:
            # #    Is there any point in this box with constraint satisfied and f < B - eps?
            # improve_formula = And(box_constraints, f_expr < B - eps)
            # m_improve = CheckSatisfiability(improve_formula, delta_dreal)
            # if m_improve is None:
            #     continue
            # else:
            #     # extract a point from m_improve, and update B
            #     intervals = [m_improve[xi] for xi in xs]
            #     mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            #     f_at_mids = rational_value_numeric(poly_num, poly_den, mids)
            #     if f_at_mids < B:
            #         B = f_at_mids
            #         # print("Updated global lower bound B =", B)

            # 2. If the box is already small enough, pick the center value
            if box.max_side_length <= min_box_size:
                f_at_mids = eval_rational(obj, box.center)
                if f_at_mids < lower_bound:
                    lower_bound = f_at_mids
                continue

            b1, b2 = split_on_longest(box)
            queue.append(b1)
            queue.append(b2)

        return lower_bound
