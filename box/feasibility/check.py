import random
from typing_extensions import Self

import numpy as np
from dreal import And, CheckSatisfiability, Expression, Formula, Not, Variable

from box.box_type import BoxN
from .type import FeasibleBox


# TODO: update with other feasibility checking methods
# perform all feasibility checking heuristics in order of execution cost
class CompleteFeasible(FeasibleBox):
    random_checks: int

    def __init__(self: Self, random_checks: int = 50):
        self.random_checks = random_checks

    def _run(
        self: Self, box: BoxN, vars: list[Variable], constr: Formula, delta: float
    ) -> bool:
        return FullFeasible()(box, vars, constr, delta)


class FullFeasible(FeasibleBox):
    def _run(
        self: Self,
        box: BoxN,
        vars: list[Variable],
        constr: Formula,
        delta: float,
    ) -> bool:
        return _full_check(box.build_constraints(vars), constr, delta)


# FIXME: type the following iface usages


class RandomFeasible(FeasibleBox):
    random_checks: int

    def __init__(self: Self, random_checks: int = 100):
        self.random_checks = random_checks

    def _run(self, box, vars, constr, delta) -> bool:
        return _random_check(box, vars, constr, self.random_checks)


class GridFeasible(FeasibleBox):
    def _run(self, box, vars, constr, delta) -> bool:
        return _grid_check(box, vars, constr, delta)


class IntervalFeasible(FeasibleBox):
    def _run(self, box, vars, constr, delta) -> bool: ...


# NOTE: dReal already implements some interval arithmetic in the solver
# TODO: learn how that works to implement better heuristics!
# TODO: how much time does this function actually take?
# check if all points in the box satisfy the given constraint
def _full_check(box_constr, constr: Formula, delta: float) -> bool:
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
def _random_check(box: BoxN, vars: list[Variable], constr: Formula, checks=100) -> bool:
    for _ in range(checks):
        pts = [random.uniform(box.min[i], box.max[i]) for i in range(len(vars))]
        # convert points to dreal Expression constants; map lazy iterator consumed on zip use
        if not constr.Substitute(dict(zip(vars, map(Expression, pts)))).Evaluate():
            return False  # found a violating point

    return True


# TODO: make this function work
# interval feasibility check, true -> feasible, false -> infeasible
def _interval_check(box: BoxN, constr: Formula) -> bool:
    # print(constr.Evaluate(box))
    # interval_eval = interval_check(box, constr)
    # if interval_eval is True:
    #     return True  # entire box is definitely feasible
    # if interval_eval is False:
    #     return False  # entire box is definitely infeasible
    # grid_check(box, constr, vars, delta)
    return True


# TODO: see if this function works
# test if entire box satisfies given constraint formula phi using a grid of test points
def _grid_check(
    box: BoxN,
    vars: list[Variable],
    constr: Formula,
    delta: float,  # 0.1
) -> bool:
    # sample density based on epsilon and box size:
    #   n_i = ceil( (max_i - min_i) / epsilon ) + 1
    samples_per_dim = [
        max(2, int(np.ceil((box.max[i] - box.min[i]) / delta)) + 1)
        for i in range(box.dim)
    ]

    # build grids
    grids = [
        np.linspace(box.min[i], box.max[i], n) for i, n in enumerate(samples_per_dim)
    ]

    # iterate over Cartesian product of grid indices
    for idx in np.ndindex(*samples_per_dim):
        pts = [float(grids[i][idx[i]]) for i in range(box.dim)]
        if not constr.Substitute(dict(zip(vars, map(Expression, pts)))).Evaluate():
            return False

    return True
