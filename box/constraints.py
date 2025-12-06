from typing import List

from dreal import And, Formula, Variable

from .type import BoxN, Point


# build dReal constraints for box
def build_constraints(box: BoxN, vars: List[Variable]) -> Formula:
    if len(vars) != box.dim:
        raise ValueError(f"Expected {box.dim} variables, got {len(vars)}")

    constraints = [
        constr
        for var, lo, hi in zip(vars, box.min, box.max)
        for constr in (lo <= var, var <= hi)
    ]

    return And(*constraints) if constraints else Formula.TRUE()


# build initial axis-aligned n dimensional box
def build_basic_box(min, max: float, dim: int) -> BoxN:
    return BoxN(Point([min] * dim), Point([max] * dim))
