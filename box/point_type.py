from dataclasses import dataclass, replace
from typing import Iterable, Iterator, Tuple

from typing_extensions import Self


# Point = NewType("Point", Tuple[float, ...])
@dataclass(frozen=True, init=False)
class Point:
    coords: Tuple[float, ...]

    def __init__(self: Self, coords: Iterable[float]):
        object.__setattr__(self, "coords", tuple(coords))

    def with_value(self: Self, i: int, v: float) -> Self:
        assert i < len(self), f"Point index {i} out of range {len(self)}."
        new = self.coords[:i] + (v,) + self.coords[i + 1 :]
        return replace(self, coords=Point(new))

    def __getitem__(self: Self, i: int) -> float:
        return self.coords[i]

    def __iter__(self: Self) -> Iterator[float]:
        return iter(self.coords)

    def __len__(self: Self) -> int:
        return len(self.coords)
