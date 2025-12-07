from abc import ABC, abstractmethod
from typing import List, Optional, Tuple

from dreal import Formula, Variable
from returns.result import Result
from typing_extensions import Self

from box import BoxN, BoxSplit, SplitLongestSide
from objective import Rational, ObjectiveBounds, AffineBounds

from .log import LogEntry


# NOTE: using abstract base class instead of protocol for explicit typing enforcement
class Algorithm(ABC):
    def _default_variables(self: Self, n: int) -> List[Variable]:
        return [Variable(f"x{i}") for i in range(n)]

    def _default_constraint(self: Self) -> Formula:
        return Formula.TRUE()

    def _default_splitter(self: Self) -> BoxSplit:
        return SplitLongestSide()

    def _default_bounder(self: Self) -> ObjectiveBounds:
        return AffineBounds()

    def __call__(
        self: Self,
        dim: int,  # explicity keep track of problem dimensionality
        init_box: BoxN,
        obj: Rational,
        vars: Optional[List[Variable]] = None,
        constr: Optional[Formula] = None,
        splitter: Optional[BoxSplit] = None,
        bounder: Optional[ObjectiveBounds] = None,
        min_box_size: float = 0.1,
        delta: float = 1e-3,
        err: float = 1e-4,
    ) -> Result[Tuple[float, List[LogEntry]], str]:
        # TODO: complete dimensionality checks
        # we want to assert that the following have the same dimensionality:
        # 1. max(# diff vars in obj fn)
        # 2. initial box dimension
        # 3. max(# diff vars in constraint formula)
        # 4. # vars in variable list

        # NOTE: would like to push these checks to type safety level,
        # but this will have to work for now
        assert init_box.dim == dim

        if vars:
            assert len(vars) == dim

        # initialize defaults for vars and constraints
        # dReal variables x0, x1, ..., x_{n-1}, w/ n = box.dim
        vars = vars or self._default_variables(dim)
        constr = constr or self._default_constraint()
        splitter = splitter or self._default_splitter()
        bounder = bounder or self._default_bounder()

        return self._run(
            dim,
            init_box,
            obj,
            vars,
            constr,
            splitter,
            bounder,
            min_box_size,
            delta,
            err,
        )

    @abstractmethod
    def _run(
        self: Self,
        dim: int,  # explicity keep track of problem dimensionality
        init_box: BoxN,
        obj: Rational,
        vars: List[Variable],  # algo impl may assume concrete vars, constr
        constr: Formula,
        splitter: BoxSplit,
        bounder: ObjectiveBounds,
        min_box_size: float,
        delta: float,
        err: float,
    ) -> Result[Tuple[float, List[LogEntry]], str]: ...
