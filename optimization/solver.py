import time
from collections import deque
from dataclasses import dataclass

from bounds import BoundingStrategy
from box import constraints, fully_feasible, max_side, split

from dreal import *


@dataclass
class OptimizationResult:
    optimal_value: float
    solve_time: float
    iterations: int
    method: str


class BranchAndBoundSolver:
    def __init__(self, bounding_strategy: BoundingStrategy):
        self.bounding_strategy = bounding_strategy
        self.iterations = 0

    def solve(
        self,
        problem: "Problem",
        delta_dreal: float = 1e-3,
        min_box_size: float = 0.1,
        eps: float = 1e-4,
    ) -> OptimizationResult:
        """
        Solve optimization problem using branch-and-bound.

        Args:
            problem: Problem instance with objective, constraints, and initial box
            delta_dreal: δ parameter for dReal precision
            min_box_size: Stop splitting when max box side < this
            eps: Pruning margin for bound comparison

        Returns:
            OptimizationResult with optimal value and metadata
        """
        start_time = time.time()
        self.iterations = 0

        # Setup dReal variables
        variables = [Variable(f"x{i}") for i in range(problem.dimension)]
        objective_expr = problem.objective(*variables)
        constraint_expr = problem.constraint(*variables)

        # Initial feasibility check
        initial_constraints = And(
            constraints(problem.initial_box, variables), constraint_expr
        )

        model = CheckSatisfiability(initial_constraints, delta_dreal)
        if model is None:
            raise ValueError("No feasible point in initial region")

        # Initialize with feasible point
        intervals = [model[var] for var in variables]
        mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
        current_best = problem.evaluate_objective(mids)

        # Main branch-and-bound loop
        queue = deque([problem.initial_box])

        while queue:
            box = queue.pop()
            self.iterations += 1

            # Strategy-specific pruning
            if self.bounding_strategy.can_prune(box, current_best, eps):
                continue

            # Standard dReal pruning check
            box_constraints = constraints(box, variables)
            improve_formula = And(
                box_constraints, constraint_expr, objective_expr < current_best - eps
            )

            m_improve = CheckSatisfiability(improve_formula, delta_dreal)
            if m_improve is None:
                continue

            # Update best solution if found better point
            intervals = [m_improve[var] for var in variables]
            mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
            candidate_value = problem.evaluate_objective(mids)
            if candidate_value < current_best:
                current_best = candidate_value

            # Check termination conditions
            if max_side(box) <= min_box_size:
                # Get local bound for small box
                local_bound = self.bounding_strategy.get_local_bound(
                    box, variables, objective_expr, constraint_expr, delta_dreal
                )
                if local_bound is not None and local_bound < current_best:
                    current_best = local_bound
                continue

            # Check if box is fully feasible
            if fully_feasible(box, constraint_expr, variables, delta_dreal):
                local_bound = self.bounding_strategy.get_local_bound(
                    box, variables, objective_expr, constraint_expr, delta_dreal
                )
                if local_bound is not None and local_bound < current_best:
                    current_best = local_bound
                continue

            # Split and continue
            b1, b2 = split(box)
            queue.append(b1)
            queue.append(b2)

        solve_time = time.time() - start_time

        return OptimizationResult(
            optimal_value=current_best,
            solve_time=solve_time,
            iterations=self.iterations,
            method=f"BranchAndBound-{type(self.bounding_strategy).__name__}",
        )


# baseline solver using dReal directly to minimize
def baseline_solve(problem: "Problem", delta_dreal: float = 1e-3) -> OptimizationResult:
    start_time = time.time()

    variables = [Variable(f"x{i}") for i in range(problem.dimension)]
    objective_expr = problem.objective(*variables)
    constraint_expr = problem.constraint(*variables)

    region_constraints = And(
        constraints(problem.initial_box, variables), constraint_expr
    )

    # Feasibility check
    model = CheckSatisfiability(region_constraints, delta_dreal)
    if model is None:
        raise ValueError("No feasible point in initial region")

    # Direct optimization
    sol_box = Minimize(objective_expr, region_constraints, delta_dreal)

    intervals = [sol_box[var] for var in variables]
    mids = [0.5 * (iv.lb() + iv.ub()) for iv in intervals]
    optimal_value = problem.evaluate_objective(mids)

    solve_time = time.time() - start_time

    return OptimizationResult(
        optimal_value=optimal_value,
        solve_time=solve_time,
        iterations=1,
        method="Baseline-dReal",
    )


def optimize(
    problem: "Problem",
    method: str = "basic",
    delta_dreal: float = 1e-3,
    min_box_size: float = 0.1,
    eps: float = 1e-4,
) -> OptimizationResult:
    """
    Args:
        problem: Problem instance
        method: "basic", "bernstein", or "baseline"
        delta_dreal: δ parameter for dReal
        min_box_size: Minimum box size for subdivision
        eps: Pruning epsilon

    Returns:
        OptimizationResult
    """
    if method == "baseline":
        return baseline_solve(problem, delta_dreal)

    # Import here to avoid circular imports
    from .bounds import BasicBounds, BernsteinBounds

    if method == "basic":
        strategy = BasicBounds()
    elif method == "bernstein":
        strategy = BernsteinBounds()
    else:
        raise ValueError(f"Unknown method: {method}")

    solver = BranchAndBoundSolver(strategy)
    return solver.solve(problem, delta_dreal, min_box_size, eps)
