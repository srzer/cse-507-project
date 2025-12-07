from abc import ABC, abstractmethod

from dreal import Formula, Variable
from typing_extensions import Self

from box.box_type import BoxN

# TODO: extract delta float to NewType Delta


class FeasibleBox(ABC):
    def __call__(
        self: Self,
        box: BoxN,
        vars: list[Variable],
        constr: Formula,
        delta: float,
    ) -> bool:
        # eventually using box.build_constraint(vars) :: Formula
        # but don't pass in just formula to allow accessing box data
        # assert both formula deal with same vars?
        return self._run(box, vars, constr, delta)

    @abstractmethod
    def _run(
        self: Self,
        box: BoxN,
        vars: list[Variable],
        constr: Formula,
        delta: float,
    ) -> bool: ...
