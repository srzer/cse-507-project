=== dReal Algorithm Overview
The dReal solver implements this relaxation to obtain a $delta$-complete algorithm for determining satisfiability:
If there is a solution to system of inequalities in the given bounded domain up to a tolerance of $delta$,
then dReal will find it; if there is no solution, dReal returns UNSAT.
Notice that if there is no solution up to a tolerance of $delta$, then there is also no solution to the exact original problem.
The complexity class is PSPACE, but in practice dReal performs very well.

dReal uses a branch-and-prune approach:
Given the initial box constraint, it will perform interval arithmetic to obtain bounds on the constraint functions to determine whether the constraints are potentially feasible.
If not, it returns UNSAT. Otherwise, it subdivides the box into smaller boxes and repeats the procedure on each sub-box:
If a sub-box is determined to be infeasible based on the interval arithmetic bounds, it is not explored further, it is ‘pruned’.
Otherwise, the sub-box is again subdivided into smaller boxes. Interval arithmetic gives imperfect bounds on a function,
but subdividing to smaller boxes with tighter interval constraints usually results in more precise bounds, which is the advantage of subdividing.
dReal will continue to either prune or subdivide until the box size is so small that the constraint functions vary by less than $delta$ within the small box,
at which point the function value is determined on that box within $delta$ precision.

