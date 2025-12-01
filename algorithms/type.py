from abc import ABC, abstractmethod
from typing import Optional, List

from dreal import Formula, Variable

from box import BoxN
from poly import Rational


# NOTE: using abstract base class instead of protocol for explicit typing enforcement
class Algorithm(ABC):
    def _default_variables(self, n: int) -> List[Variable]:
        return [Variable(f"x{i}") for i in range(n)]

    def _default_constraint(self) -> Formula:
        return Formula.TRUE()

    def __call__(
        self,
        dim: int,  # explicity keep track of problem dimensionality
        init_box: BoxN,
        obj: Rational,
        vars: Optional[List[Variable]] = None,
        constr: Optional[Formula] = None,
        min_box_size: float = 0.1,
        delta: float = 1e-3,
        err: float = 1e-4,
    ) -> Optional[float]:
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

        self._run(dim, init_box, obj, vars, constr, min_box_size, delta, err)

    @abstractmethod
    def _run(
        self,
        dim: int,  # explicity keep track of problem dimensionality
        init_box: BoxN,
        obj: Rational,
        vars: List[Variable],  # algo impl may assume concrete vars, constr
        constr: Formula,
        min_box_size: float,
        delta: float,
        err: float,
    ) -> Optional[float]: ...
