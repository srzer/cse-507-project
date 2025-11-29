"""
Unified optimization framework for constrained rational function minimization.

This package provides a unified branch-and-bound solver that works for any dimension
with pluggable bounding strategies including basic dReal bounds and Bernstein polynomial bounds.
"""

from .solver import optimize, BranchAndBoundSolver
from .bounds import BasicBounds, BernsteinBounds
from .problems import Problem

__all__ = ['optimize', 'BranchAndBoundSolver', 'BasicBounds', 'BernsteinBounds', 'Problem']