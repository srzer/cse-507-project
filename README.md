# cse-507-project

For deal installation, please refer to https://github.com/dreal/dreal4.

## Current implementation of n-dim Brach-and-prune algo with berstein method

- The input rational function is represented as two lists of terms, each term is [coefficient, [exponent_vector]]. For example, $f(x,y) = (2x^2y + 3y^2) / (x^2 + y)$, then numerator_terms = [ [2, [2,1]], [3, [0,2]] ], and denominator_terms = [ [1, [2,0]], [1, [0,1]] ].
- We have an initial box, which is assumed to cover all feasible areas and contain no pole. We also have a minimal box size, which is set default as same as $\delta$ in dReal (e.g. 1e-3) by now.
- For each box, the pipeline is:
    - We solve (Bernstein method) a rough lower bound on this box. If it's larger than current lower bound, skip.
    - We then solve (dReal) a solution in the box satisfying the constraint and improving the lower bound at least by $\epsilon$, if unsat, skip; if sat, update the lower bound.
    - If the box is already smaller than the minimal box size, skip.
    - We then determine (dReal) whether the box is fully feasible. If so, we solve (dReal) the accurate lower bound with only box constraint; otherwise, we split the box (currently splitting the longest edge).

To do:
- The splitting can be guided using heuristic;
- Find how to replace some expensive dReal solver using Bernstein method;
- Write more test-case to show our method's advantage.