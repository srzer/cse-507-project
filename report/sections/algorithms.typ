== Our Branch & Bound Algorithm
Our general algorithmic idea is a branch-and-prune approach, as inspired by dReal:
Given the initial box constraint, we will first perform interval arithmetic to obtain bounds on the constraint functions to determine whether the constraints are potentially feasible.
If not, it returns UNSAT. Otherwise, it subdivides the box into smaller boxes and repeats the procedure on each sub-box:
If a sub-box is determined to be infeasible based on the interval arithmetic bounds, it is not explored further, it is pruned.
Otherwise, the sub-box is again subdivided into smaller boxes. Interval arithmetic gives imperfect bounds on a function,
but subdividing to smaller boxes with tighter interval constraints usually results in more precise bounds, which is the advantage of subdividing.
We will continue to either prune or subdivide until the box size is so small that the constraint functions vary by less than $delta$ within the small box,
at which point the function value is determined on that box within $delta$ precision. Now we present how we improved the vanilla idea through implementing new techniques as below.

=== Objectives
Our main focus will be rational objective optimization with constraints.
The rational objective function $f(x) = p(x)\/q(x)$ is represented by two lists of terms: the numerator terms $[ [c_i, e_i], dots ]$,
and the denominator terms $[ [d_j, g_j], dots ]$. For example, if the objective is $f(x, y) = (2x^2 y + 3y^2) \/ (x^2 + y)$,
then we have the numerator terms as $[ [2, [2,1]], [3, [0,2]] ]$, and the denominator terms as $[ [1, [2,0]], [1, [0,1]] ]$.

=== Step 0. Initialization
The algorithm begins with an initial search box that is assumed to contain all feasible solutions and does not include any poles of the rational function.
A minimal box size is also defined, usually equal to the numerical tolerance used by dReal (for example, `1e-3`).
The global lower bound is initialized, and the initial box is pushed into a queue of boxes to process.

=== Step 1. Basic pruning
For each box, the algorithm first applies a pruning step using the affine method (we've also tried Bernstein method, which turns out to be efficient in $n=3$ case,
but inefficient for larger dimensions) to estimate a rough lower bound of the objective function on that box.
If this rough bound is already worse than the current global lower bound, the box is discarded without further work.
Otherwise, the algorithm asks dReal to search for a feasible point inside the box that improves the current lower bound by at least a small threshold.
If dReal reports that no such point exists, the box is discarded. If a better point is found, the global lower bound is updated.

=== Step 2. Handling full feasibility
If the box is already smaller than the minimal box size, it is skipped. Otherwise, the algorithm uses dReal to check whether the entire box is feasible.
If it is fully feasible, then we provide two approaches to handle this.
The approach 1 is to use dreal’s Minimize to directly obtain a lower bound,
note that the box is smaller than the whole original box, so it is still efficient;
the approach 2 is to invoke an iterative splitting process to obtain a lower bound on this box using affine method.

=== Step 3. Splitting
If the box is only partially feasible, the algorithm simply splits it into two smaller boxes and continues processing.
Beyond simply splitting the longest edge into half, which is the default implementation, we implemented a gradient-guided splitting heuristic,
which calculates the gradient direction of the objective to optimize, and then selects the edge which best aligns the direction of the gradient.
