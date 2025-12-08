There are two key elements to introduce: Interval Arithmetic and the dReal SMT solver. @gao2013dreal
Given a collection of variables $x_j$ whose values are each constrained to given intervals $I_j = [a_j, b_j]$, the problem of interval arithmetic is determining the most precise interval to which some arithmetic combination of these variables are constrained. For instance, given $x$ and $y$, both constrained to the interval $[0, 1]$, we can say that $x + y$ is constrained to the interval $[0, 2]$ and $x - y$ is constrained to the interval $[-1, 1]$.
Complications emerge when there are dependencies between our variables;
for instance, supposing we have the relation $x=y$, then $x - y$ is actually constrained to the interval $[0, 0]$.
We can restate this problem as follows: given a polynomial on a box $f: I_1 times I_2 times dots.c times I_m -> RR$,
determine upper and lower bounds on the value of $f$ subject to constraints $g_i (x_1, x_2, dots, x_m) >= 0$.

There are a number of approaches to interval arithmetic.
The goal is to find bounds that are as tight as possible, as efficiently as possible.
In general, different approaches have different strengths and weaknesses, and no single method dominates all others across all test cases.
Various types of branch-and-prune methods, such as those described below, are used, along with different heuristics for search and subdivision.
Some methods are explicitly designed for polynomials, while others can generalize to transcendental functions.
#footnote[Some algorithms handle transcendental functions by simply using Taylor approximations.]
Given any set of inequalities involving arithmetic combinations and compositions of polynomials and transcendental functions in several variables (a very broad scope),
we are interested in determining whether there are some real assignments of the input variables that will simultaneously satisfy every given inequality.
Note equalities can be given as a pair of inequalities, and so are also included.
An example problem would be: Is there some real $x, y$ satisfying $sin(x) + y = e^(x y) and tan(x^2 + y^2) >= log(x - y)$?
These are difficult problems;
in fact, generically they are undecidable.
This means we cannot simply feed these problems into an SMT solver using nonlinear real arithmetic with transcendental functions.

However, the problem becomes decidable if we take two actions:
First, we want to restrict our search to a compact domain, say a box, by putting bounds on every variable individually.
Second, we can introduce an error tolerance instead of looking for an exact solution.
That is, we can choose any $delta > 0$, and then for every constraint $f(x) >= 0$, we weaken this to $f(x) >= -delta$,
and an equality $f(x) = 0$ will become $|f(x)| <= delta$.
