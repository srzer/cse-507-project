from abc import ABC, abstractmethod
from typing import Optional, Tuple
from typing_extensions import Self

from objective import Rational, eval_gradient

from box import BoxN
from objective.polynomial_type import unit

# NOTE: rewrite these as Box methods using mixins?


# TODO: create splitting heuristsics class to choose
# and plug-in different splits
class BoxSplit(ABC):
    def _default_objective(self: Self, n: int) -> Rational:
        return Rational(unit(n), unit(n))

    def __call__(self: Self, box: BoxN, obj: Optional[Rational]) -> Tuple[BoxN, BoxN]:
        obj = obj or self._default_objective(box.dim)
        return self._run(box, obj)

    @abstractmethod
    def _run(self: Self, box: BoxN, obj: Rational) -> Tuple[BoxN, BoxN]: ...


class SplitLongestSide(BoxSplit):
    def _run(self: Self, box: BoxN, obj: Rational) -> Tuple[BoxN, BoxN]:
        return _split_on_longest(box)


class SplitGradient(BoxSplit):
    def _run(self: Self, box: BoxN, obj: Rational) -> Tuple[BoxN, BoxN]:
        return _split_by_gradient(obj, box)


# TODO: improve function doc
# split in direction of greatest function slope?
def _split_by_gradient(f: Rational, box: BoxN) -> Tuple[BoxN, BoxN]:
    grad = eval_gradient(f, box.center)
    nonzero_gradient = any(abs(x) > 0.0 for x in grad)
    # use box interval lengths as default measure, which is the same as split on longest
    box_side_measures = list(map(abs, grad)) if nonzero_gradient else box.lengths
    max_side_idx = max(range(box.dim), key=lambda i: box_side_measures[i])
    mid = box.center[max_side_idx]
    return (
        BoxN(box.min, box.max.with_value(max_side_idx, mid)),
        BoxN(box.min.with_value(max_side_idx, mid), box.max),
    )


# split on longest dimension
def _split_on_longest(box: BoxN) -> Tuple[BoxN, BoxN]:
    idx = box._max_side_idx()
    mid = box.center[idx]

    return (
        BoxN(box.min, box.max.with_value(idx, mid)),
        BoxN(box.min.with_value(idx, mid), box.max),
    )
