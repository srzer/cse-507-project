from abc import ABC, abstractmethod
from typing_extensions import Self

from box import BoxN
from interval import Bounds

from objective.rational_type import Rational

# TODO: implement these as defaults in the abstract class
ZERO_TOLERANCE = 1e-15  # tolerance for filtering zero coefficients
DENOM_MIN_THRESHOLD = 1e-8  # min demoninator to prevent div-by-zero


# NOTE: this class has objective as first arg
# instead of BoxSplit, which uses box as first,
# but both classes take a pair of obj and box.
class ObjectiveBounds(ABC):
    def __call__(
        self: Self,
        obj: Rational,
        box: BoxN,
    ) -> Bounds:
        if obj.n_vars != box.dim:
            raise ValueError(
                f"Objective function of {obj.n_vars} mismatch with box of {box.dim} dimensions."
            )

        return self._run(obj, box)

    @abstractmethod
    def _run(self: Self, obj: Rational, box: BoxN) -> Bounds: ...
