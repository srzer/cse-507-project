import random

import numpy as np
from dreal import And, CheckSatisfiability, Expression, Formula, Not, Variable

from .type import BoxN


# implements all feasibility checking heuristics in order of execution cost
def feasible(
    box: BoxN,
    constr: Formula,
    box_constr: Formula,
    delta: float,
    random_checks: int = 50,
) -> bool:
    # interval_eval = interval_check(box, constr)
    # if interval_eval is True:
    #     return True  # entire box is definitely feasible
    # if interval_eval is False:
    #     return False  # entire box is definitely infeasible
    # grid_check(box, constr, vars, delta)
    return full_check(constr, box_constr, delta)


# check if all points in the box satisfy the given constraint
def full_check(constr: Formula, box_constr: Formula, delta: float) -> bool:
    """
    The whole box is feasible iff
    there is NO point in the box violating the constraint.

    box is fully feasible  ⇔  (box_constraints ∧ ¬constraint) is UNSAT.
    """
    formula_violate: Formula = And(box_constr, Not(constr))
    model = CheckSatisfiability(formula_violate, delta)
    # UNSAT ⇒ no violating point ⇒ box fully feasible
    # result = constraint.Evaluate(box)
    return model is None


# TODO: implement this function properly
# random feasibility check with n points
def random_check(box: BoxN, constr: Formula, vars: list[Variable], n=100) -> bool:
    for _ in range(n):
        pts = [random.uniform(box.min[i], box.max[i]) for i in range(len(vars))]
        # convert points to dreal Expression constants; map lazy iterator consumed on zip use
        if not constr.Substitute(dict(zip(vars, map(Expression, pts)))).Evaluate():
            return False  # found a violating point

    return True


# TODO: make this function work
# interval feasibility check, true -> feasible, false -> infeasible
def interval_check(box: BoxN, constr: Formula):
    # print(constr.Evaluate(box))
    return True


# TODO: see if this function works
# test if entire box satisfies given constraint formula phi using a grid of test points
def grid_check(
    box: BoxN,
    constr: Formula,
    vars: list[Variable],  # FIXME: why do we have this?
    epsilon: float = 0.1,
) -> bool:
    # sample density based on epsilon and box size:
    #   n_i = ceil( (max_i - min_i) / epsilon ) + 1
    samples_per_dim = [
        max(2, int(np.ceil((box.max[i] - box.min[i]) / epsilon)) + 1)
        for i in range(box.dim)
    ]

    # build grids
    grids = [
        np.linspace(box.min[i], box.max[i], n)
        for i, n in enumerate(samples_per_dim)
        # weird generator type, want to cast to float
    ]

    # iterate over Cartesian product of grid indices
    for idx in np.ndindex(*samples_per_dim):
        pts = [float(grids[i][idx[i]]) for i in range(box.dim)]
        if not constr.Substitute(dict(zip(vars, map(Expression, pts)))).Evaluate():
            return False

    return True
