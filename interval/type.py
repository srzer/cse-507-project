from dataclasses import dataclass
from typing import Iterator
from typing_extensions import Self
# from box import Point
# from dreal import Box, Interval, Variable


@dataclass(frozen=True)
# numeric bounds for a single dimension [low, high]
class Bounds:
    lo: float
    hi: float

    def __post_init__(self: Self):
        if self.lo > self.hi:
            raise ValueError(f"Invalid bounds: {self.lo} > {self.hi}")

    @property
    def width(self: Self) -> float:
        return self.hi - self.lo

    @property
    def mid(self: Self) -> float:
        return (self.lo + self.hi) / 2

    def __iter__(self: Self) -> Iterator[float]:
        return iter([self.lo, self.hi])


# def from_model(model: Box, vars: list[Variable]) -> list[Interval]:
#     return [model[x_i] for x_i in vars]


# get middle of interval list
# def mids(ivs: list[Interval]) -> Point:
#     return Point(tuple([iv.mid() for iv in ivs]))
