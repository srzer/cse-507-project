from dataclasses import dataclass
from typing import List

from dreal import Formula, Variable

from box import BoxN
from objective import Rational


@dataclass(frozen=True)
# dataclass to represent a test problem for algos
class Problem:
    name: str
    dim: int
    objective: Rational
    constraints: Formula
    initial_box: BoxN
    variables: List[Variable]
