from collections import deque
from typing import Deque, List
from typing_extensions import Self

from dreal import Formula, Variable

from box import BoxN, BoxSplit, Point
from objective import ObjectiveBounds, Rational, eval_rational

from .either import Either, Left, Right
from .type import Algorithm

# FIXME: why do we have separate algorithm implementations?
# shouldn't we just have all the heurstic pieces (or other components)
# just be passables to plug and play?
# ... i suppose we can do this for now


class FeasibleMinBranchAndBound(Algorithm):
    initial_lower_bound: float

    def __init__(self: Self, initial_lower_bound: float = float("inf")):
        self.initial_lower_bound = initial_lower_bound

    def _run(
        self: Self,
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
        # fn_expr = eval_symbolic(obj, vars)
        lower_bound = self.initial_lower_bound

        queue: Deque[BoxN] = deque()
        queue.append(init_box)

        while queue:
            box: BoxN = queue.pop()

            bounded_min, _ = bounder(obj, box)

            if bounded_min >= lower_bound - err:
                continue

            # NOTE: corners scale exponentially with dimension
            # TODO: implement heuristic to prevent over-execution
            for corner in generate_corners(box):
                obj_at_corner = eval_rational(obj, corner)
                if obj_at_corner < lower_bound:
                    lower_bound = obj_at_corner

            # 2. If the box is already small enough, pick the center value
            if box.max_side_length <= min_box_size:
                obj_at_center = eval_rational(obj, box.center)
                if obj_at_center < lower_bound:
                    lower_bound = obj_at_center
                continue

            b1, b2 = splitter(box, obj)
            queue.append(b1)
            queue.append(b2)

        return Right(lower_bound)


# generate all 2^n corner points of a box in n dimensions
def generate_corners(box: BoxN) -> list[Point]:
    corners = []

    def _generate_recursive(idx: int, current: list[float]) -> None:
        if idx == box.dim:
            corners.append(Point(tuple(current)))
            return

        # try both low and high values for this dimension
        current.append(box.min[idx])
        _generate_recursive(idx + 1, current)
        current.pop()

        current.append(box.max[idx])
        _generate_recursive(idx + 1, current)
        current.pop()

    _generate_recursive(0, [])
    return corners
