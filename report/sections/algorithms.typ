== Our Branch & Bound Algorithm
Our algorithm follows a branch-and-prune strategy inspired by dReal.

Starting from the initial box constraint, we use interval arithmetic to bound the constraint functions
and determine whether the box may contain feasible points. If the box is infeasible, we return UNSAT.
Otherwise, we subdivide the box and apply the same procedure recursively.

Any sub-box that is proven infeasible by interval arithmetic is pruned.
Sub-boxes that remain feasible are further subdivided.
Although interval arithmetic provides loose bounds, working on smaller boxes tightens these bounds
and makes the feasibility test more accurate.

This process continues until each box becomes sufficiently small so that the constraint functions vary
by less than $δ$ inside the box. At that point, the function value on the box is determined within $δ$ precision.
Below, we describe how we enhance this basic framework with several additional techniques.


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

=== Step 3. Splitting heuristic
If the box is only partially feasible, the algorithm simply splits it into two smaller boxes and continues processing.
Beyond simply splitting the longest edge into half, which is the default implementation, we implemented a gradient-guided splitting heuristic,
which calculates the gradient direction of the objective to optimize, and then selects the edge which best aligns gradient direction.
