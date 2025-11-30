from dataclasses import dataclass
from box import Point
from dreal import Box, Interval, Variable

# @dataclass(frozen=True)
# class Interval:
#     min: float
#     max: float


# def from_model(model: Box, vars: list[Variable]) -> list[Interval]:
#     return [model[x_i] for x_i in vars]


# get middle of interval list
# def mids(ivs: list[Interval]) -> Point:
#     return Point(tuple([iv.mid() for iv in ivs]))
