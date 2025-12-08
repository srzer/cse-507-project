#import "@preview/elsearticle:1.1.0": *

== Summary
Performance differences may derive from the structure of the search space.
Our method works well on some problems, due to our specialization on rational functions.
On other problems, dReal stood out compared with all other surveyed methods.

On Pole Avoidance in Figure @fig:problem_runtime_pole_avoidance:
The objective function is $1 \/ (x + y + z - 2.5)$ on a simple box.
This function is monotonic and smooth over the entire feasible region.
The gradient provides a very strong and consistent signal that our branch-and-bound algorithm can use effectively to quickly prune large parts of the space and converge on the minimum without much wasted effort.

Our method is significantly faster than the baseline (normalized runtime is fast, around 0.01-0.02x the baseline time)
and numeric optimizations and finds a slightly better lower bound.

#figure(
  image("../renders/full_comparison_Pole_Avoidance.svg", width: 100%),
  caption: [Pole Avoidance runtime.
    All of our methods shine. Due for further exploration.
  ],
) <fig:problem_runtime_pole_avoidance>

On Positive Islands in Figure @fig:problem_runtime_dreal_outlier:
Positive Islands The objective function is a simple polynomial $x^2 + 1$,
but the feasible region consists of two small, disconnected spheres in the middle of a very large search box.
Our algorithm likely struggles here because it has to spend the vast majority of its time subdividing a large, empty, infeasible space just to find the two tiny islands where solutions exist.
The dReal baseline method might be more efficient at handling these kinds of disjoint feasible regions,
which would explain why our method is comparatively slower on this specific problem.

#subfigure(
  figure(
    image("../renders/full_comparison_Positive_Islands.svg", width: 100%),
    caption: [Positive islands runtime performance.],
  ),
  <fig:problem_runtime_positive_islands>,
  figure(
    image("../renders/full_comparison_Singularity_Edge.svg", width: 100%),
    caption: [Singularity edge runtime performance],
  ),
  <fig:problem_runtime_singularity_edge>,
  columns: (1fr, 1fr),
  caption: [(a) dReal defeats all other methods with marginal error.
    SHGO and Dual Annealing get the wrong answer.
    (b) Our method is able to match standard numerical optimization performance.

  ],
  label: <fig:problem_runtime_dreal_outlier>,
)

In @fig:problem_runtime_positive_islands, our method is much slower than the dReal baseline (normalized runtime is around 6-7x slower the baseline time) while finding a similar bound.
This is also the case for numeric optimizers.



=== Aggregate Runtime Performance

#figure(
  image("../renders/aggregate_normalized_runtime.svg", width: 90%),
  caption: [Aggregate runtime performance across testing suite.
    Splitting on longest box side yielded a faster median than baseline.],
) <fig:aggregate_runtime>

We now consider the aggregate of runtimes across the testing suite.
Our method outperformed the dReal baseline in speed, but generally found smaller lower bounds.

#figure(
  image("../renders/full_aggregate_normalized_runtime.svg", width: 90%),
  caption: [Plot of aggregate runtime performance across testing suite with standard numeric optimization methods.
    Note that the order methods differs from @fig:aggregate_runtime.
    All numerical methods have a faster median runtime.],
) <fig:aggregate_runtime_numeric>

However, runtime cannot be considered in isolation of the actual derived bound.
We see that all surveyed methods tend to find a smaller mininum than the dReal baseline, which happens to be incorrect for some of the problems.

#figure(
  image("../renders/full_aggregate_bound_difference.svg", width: 90%),
  caption: [Normalized aggregate bound difference across testing suite.
  ],
) <fig:aggregate_bound_numeric>


//TODO: would be great to get visualizations on each of the particular problem plots

== Design & Implementation Challenges

We discovered that the affine bounding heuristic can perform rapid minimization, but the heuristic alone is insufficient to uncover the actual minimum.
We consider the following Split Islands problem.
The algorithms with the affine interval bounding heuristic terminate quickly, but with the incorrect lower bound.
#figure(
  image("../renders/full_comparison_Split_Islands.svg", width: 90%),
  caption: [Runtime for Split Islands.],
) <fig:problem_runtime_split_islands>


For some problems, numerical methods outperform solver-aided ones.
The following plots demonstrate a coupling between our method and dReal with respect to runtime performance on the Rational Valley and Sanity Poly problems.

#subfigure(
  figure(
    image("../renders/full_comparison_Rational_Valley.svg", width: 100%),
    caption: [Rational valley runtime performance.],
  ),
  <fig:problem_runtime_rational_valley>,
  figure(
    image("../renders/full_comparison_Sanity_Poly.svg", width: 100%),
    caption: [Sanity poly runtime performance],
  ),
  <fig:problem_runtime_sanity_poly>,
  columns: (1fr, 1fr),
  caption: [(a) On certain problems, solvers have a difficult time.
    (b) On others, the numerical methods suffer.
    All optimizers agree on the lower bound for both problems.
  ],
  label: <fig:problem_runtime_numeric_comparison>,
)


// == Design & Implementation Challenges



== Selecting Bencmarks

We collected sample benchmarks to form our test suite from Wikipedia. #footnote[https://en.wikipedia.org/wiki/Test_functions_for_optimization]
We currently have tests for:
Sanity Poly, Sanity Rational, Rational Bowl, Rational Valley, Split Islands, Positive Islands, Singularity Edge, Pole Avoidance, Sparse Intersection, and also Himmelblau Ratio (but only for 3D).

// == Improvement Areas

// Adaptive Heuristics:
// Our results show that some heuristic combinations are fast but coarse, while others are slow but precise.
// We could define transition conditions from one splitting and bounding combo to another.
// #footnote[These conditions would be another heuristic themselves. We might wish to avoid second order heuristics.]
// For example, an adaptive solver could begin with a prospecting phase using a fast heuristic like affine bounds with longest side bisect splitting.
// This would rapidly discard large regions of the search space and establish an initial upper bound on the global minimum.
// Once the rate of improvement slows down or the total volume of the active search boxes falls below a threshold,
// the algorithm could switch to a refining phase, using a more expensive but precise heuristic like the Bernstein bounds and gradient split to zero in on the true minimum.

// Advanced Box Splitting Strategies:
// The current project uses bisecting on the longest box side which is ignorant to function behaviour,
// and splitting in the direction of greatest change in the objective function, gradient split.
// There's a rich area of research here that could yield significant performance gains.
// Additionally, it would be beneficial to compare our attempts at smart methods to merely random splitting.
// Another intelligent heuristic for non-linear functions could be to split where the function has the greatest measure,that is where it is changing the most.
// We could implement a maximal smear splitter. For each dimension of a box, it would calculate the width of the objective function's range (its “smear”), an interval.
// #footnote[So, we could apply interval arithmetic SMT tools here.]
// It would then split the box along the dimension that produces the widest range, as this is often where the most progress can be made in tightening the bounds.
// This could be particularly effective on problems like Rational Valley, where the function's behavior is complex and not well-aligned with the box's geometry.

// Ordering Hybrid SMT and Numerical Methods.
// We treat the SMT and numerical methods as separate categories for comparison, but their strengths are complementary.
// Our hybrid approach of refining the search space early with numerical methods before making calls to dReal aims to be more powerful than either technique alone.
// However, dReal does have similar baked-in optimizations itself, and is not purely an SMT tool.
// Our current approach places the SMT-based dReal at the end of the pipeline,
// yet we could reverse this by using the SMT-based branch-and-bound algorithm to find and isolate guaranteed feasible regions before the final minimum.
// The SMT solver is excellent at handling complex constraints, as seen in the Positive Islands problem.
// Once it identifies a small box that is provably feasible, it could hand that box off to a fast numerical optimizer
// like SHGO or Differential Evolution to quickly find the local minimum within that feasible region.


