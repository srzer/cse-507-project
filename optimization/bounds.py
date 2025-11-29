"""
Bounding strategies for the branch-and-bound optimization algorithm.

This module provides different strategies for computing bounds and pruning
the search space during optimization.
"""

from typing import List, Optional, Protocol
from abc import ABC, abstractmethod

from dreal import *

from box import Box, constraints


# Defines requirements for different bounding strategies
class BoundingStrategy(ABC):
    @abstractmethod
    def can_prune(self, box: Box, current_bound: float, eps: float) -> bool:
        # return True when given box can be pruned based on bounds
        ...

    @abstractmethod
    def get_local_bound(
        self,
        box: Box,
        variables: list,
        objective_expr: Formula,
        constraint_expr: Formula,
        delta: float,
    ) -> Optional[float]:
        # get local bound for given box, otherwise None if infeasible
        ...


# basic bounding only using dReal minimizing
class BasicBounds(BoundingStrategy):
    """Basic bounding strategy using only dReal minimize."""

    def can_prune(self, box: Box, current_bound: float, eps: float) -> bool:
        """Basic strategy has no additional pruning beyond what's in the main loop."""
        return False

    def get_local_bound(
        self,
        box: Box,
        variables: List[Variable],
        objective_expr: Formula,
        constraint_expr: Formula,
        delta: float,
    ) -> Optional[float]:
        """Get local bound using dReal minimize."""
        box_constraints = constraints(box, variables)
        local_constraints = And(box_constraints, constraint_expr)

        try:
            sol_box = Minimize(objective_expr, local_constraints, delta)
            if sol_box is None:
                return None

            # Get midpoint and evaluate objective
            intervals = [sol_box[var] for var in variables]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]

            # Convert symbolic objective to numerical evaluation
            # We need to substitute variables with values
            subst = {var: val for var, val in zip(variables, mids)}
            result = objective_expr.Substitute(subst).Evaluate()

            return float(result)
        except:
            return None


class BernsteinBounds:
    """Bernstein polynomial bounding strategy."""

    def __init__(self):
        try:
            from bernstein_utils import bernstein_bounds_on_box
            from poly_utils import poly_from_terms

            self.bernstein_bounds_on_box = bernstein_bounds_on_box
            self.poly_from_terms = poly_from_terms
        except ImportError:
            print(
                "Warning: Bernstein utilities not available, falling back to basic bounds"
            )
            self.bernstein_bounds_on_box = None

    def can_prune(self, box: Box, current_bound: float, eps: float) -> bool:
        """Use Bernstein polynomial bounds for pruning."""
        if self.bernstein_bounds_on_box is None:
            return False

        try:
            # This would need to be adapted based on the specific problem
            # For now, we'll disable Bernstein pruning and focus on the architecture
            return False
        except Exception:
            return False

    def get_local_bound(
        self,
        box: Box,
        variables: List[Variable],
        objective_expr: Formula,
        constraint_expr: Formula,
        delta: float,
    ) -> Optional[float]:
        """Get local bound using both Bernstein bounds and dReal minimize."""
        # Fall back to basic bounds for now
        basic_bounds = BasicBounds()
        return basic_bounds.get_local_bound(
            box, variables, objective_expr, constraint_expr, delta
        )


# class EnhancedBounds:
#     """Enhanced bounding strategy that combines multiple techniques."""

#     def __init__(self):
#         self.basic = BasicBounds()
#         self.bernstein = BernsteinBounds()

#     def can_prune(self, box: Box, current_bound: float, eps: float) -> bool:
#         """Try Bernstein pruning first, then fall back to basic."""
#         return self.bernstein.can_prune(
#             box, current_bound, eps
#         ) or self.basic.can_prune(box, current_bound, eps)

#     def get_local_bound(
#         self,
#         box: Box,
#         variables: List[Variable],
#         objective_expr: Formula,
#         constraint_expr: Formula,
#         delta: float,
#     ) -> Optional[float]:
#         """Use the best available bound."""
#         bernstein_bound = self.bernstein.get_local_bound(
#             box, variables, objective_expr, constraint_expr, delta
#         )
#         basic_bound = self.basic.get_local_bound(
#             box, variables, objective_expr, constraint_expr, delta
#         )

#         if bernstein_bound is None:
#             return basic_bound
#         elif basic_bound is None:
#             return bernstein_bound
#         else:
#             return min(bernstein_bound, basic_bound)
