If we have a way to determine satisfiability of a set of constraints, we also have a way to perform optimization:
To minimize a function $f(x)$ subject to constraints $g_i (x)$ on a given box domain,
we can repeatedly check the satisfiability of the system $f(x) <= C$ together with $g_i (x)$.
If this system is satisfiable, we know the optimal solution is at most $C$, and if it is not satisfiable then the optimal solution is at least $C$.
Therefore if we start with a known lower bound on the value of $f(x)$ and a known feasible solution (to get an upper bound on the minimum),
then using a binary search will efficiently converge to the solution of the optimization problem. Using this strategy, dReal can solve optimization problems.

We can hope to make this more efficient by performing this bounding of $f(x)$ simultaneously to our branch-and-prune optimization procedure.
Interval arithmetic on $f(x)$ will give us an initial lower bound on $f(x)$.
As we branch-and-prune in a depth-first search to find regions satisfying the constraints and eliminate regions that are infeasible,
once we find a region where the constraints are δ-satisfiable we can find an upper bound for $f(x)$ in this small region.
This gives an upper bound on our optimal solution. From this point forwards, in our branch-and-prune process we have another mechanism for pruning:
Using interval arithmetic, we can find a lower bound for $f(x)$ on a given box, and if this lower bound is greater than our current upper bound on the minimum of $f(x)$,
then we know $f(x)$ will not be optimized within this box and we can prune the entire box.
Eventually we will explore the entire feasible domain and be left with small feasible regions where $f(x)$ is smallest.
Having reduced the problem to small boxes, it is then efficient to find the minimum of $f(x)$ on each of these boxes directly using the method described above, as implemented in dReal.

// maybe including a figure here to keep things interesting
#figure(
  image("../renders/examples_2.png", width: 90%),
  caption: [From left to right: hand-written example, Singular Edge, and Rational Bowl optimization search space renders.
  ],
) <fig:sample_problem_renders>

Our goal is to see if this approach can outperform dReal’s approach to optimization.
In particular, we considered rational functions in several variables, and hoped to optimize the procedure for this specific class of functions.
In addition to the use of dReal, our approach drew heavily from two previous techniques: Bernstein polynomials and affine arithmetic.

=== Bernstein Polynomiial Form
Bernstein polynomials give an algorithm for computing exact upper and lower bounds on a given polynomial.
It is computationally expensive to compute the initial bounds,
however the method of Bernstein polynomials interfaces well with subdividing along boxes so that subsequent refined bounds are more efficient to compute.
We apply this method separately on the numerator and denominator of our objective function and then use naive interval arithmetic to bound the resulting ratio, which results in an inexact bound.

=== Affine Arithmetic
Affine arithmetic is a method for obtaining an approximate interval bound on a polynomial,
and there is both a partial/first-order affine arithmetic and a full/higher-order affine arithmetic, the latter being more precise but involving more computation.
Bernstein polynomials are very technical, but the method of affine arithmetic can be illustrated simply by example.
If our variables are constrained to the intervals $x in [2, 4]$ and $y [3,7]$, then we define auxiliary variables constrained to $[-1, 1]$,
call them $e_1$ and $e_2$, so that $x = e_1 + 3$ and $y = 2 e_2 + 5$. This normalization allows us to handle dependencies;
when we compute the interval for $x -$x we compute $(e_1 + 3) - (e_1 + 3) = 0$.
The purpose of normalizing to $[-1, 1]$ is so that products of these auxiliary variables are again variables lying in $[-1, 1]$.
For instance, when applying affine arithmetic to the function $x y$ we will get a cross term with $e_1 e_2$.
In the partial or first-order affine arithmetic, we would simply substitute in $-1$ or $1$ here when computing extreme values.
In higher-order affine arithmetic, these cross terms are replaced with new auxiliary variables, and then again dependencies between different instances of these cross terms can be tracked.

== Our Branch & Bound Algorithm
Our general algorithmic idea is a branch-and-prune approach, which is inspired by dReal:
Given the initial box constraint, it will perform interval arithmetic to obtain bounds on the constraint functions to determine whether the constraints are potentially feasible.
If not, it returns UNSAT. Otherwise, it subdivides the box into smaller boxes and repeats the procedure on each sub-box:
If a sub-box is determined to be infeasible based on the interval arithmetic bounds, it is not explored further, it is pruned.
Otherwise, the sub-box is again subdivided into smaller boxes. Interval arithmetic gives imperfect bounds on a function,
but subdividing to smaller boxes with tighter interval constraints usually results in more precise bounds, which is the advantage of subdividing.
dReal will continue to either prune or subdivide until the box size is so small that the constraint functions vary by less than $delta$ within the small box,
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
For each box, the algorithm first applies a pruning step using the affine method (we’ve also tried Bernstein method, which turns out to be efficient in $n=3$ case,
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
