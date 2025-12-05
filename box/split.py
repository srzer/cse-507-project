from typing import Tuple

from poly import Rational

from box import BoxN

from poly.type import eval_gradient

# NOTE: rewrite these as Box methods using mixins?


# TODO: create splitting heuristsics class to choose
# and plug-in different splits


# TODO: improve function doc
# split in direction of greatest function slope?
def split_box_gradient(f: Rational, box: BoxN) -> Tuple[BoxN, BoxN]:
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
def split_on_longest(box: BoxN) -> Tuple[BoxN, BoxN]:
    idx = box._max_side_idx()
    mid = box.center[idx]

    return (
        BoxN(box.min, box.max.with_value(idx, mid)),
        BoxN(box.min.with_value(idx, mid), box.max),
    )
