#import "@preview/elsearticle:1.1.0": *

== Summary
Performance differences may derive from the structure of the search space.
Our method works well on some problems, perhaps due to our specialization on ration functions.
On other problems, dReal an outlier compared with all other methods.

On Pole Avoidance in @fig:problem_runtime_pole_avoidance:
The objective function is $1 \/ (x + y + z - 2.5)$ on a simple box.
This function is monotonic and smooth over the entire feasible region.
The gradient provides a very strong and consistent signal that our branch-and-bound algorithm can use effectively to quickly prune large parts of the space and converge on the minimum without much wasted effort.


#figure(
  image("../renders/full_comparison_Pole_Avoidance.svg", width: 110%),
  caption: [Plot of runtime performance for Pole Avoidance.
    Our method is significantly faster than the baseline (normalized runtime is very low, around 0.01-0.02x the baseline time)
    and numeric optimizations and finds a slightly better (lower) bound.
  ],
) <fig:problem_runtime_pole_avoidance>

On Positive Islands in @fig:problem_runtime_dreal_outlier:
Positive Islands The objective function is a simple polynomial $x^2 + 1$,
but the feasible region consists of two small, disconnected spheres in the middle of a very large search box.
Our algorithm likely struggles here because it has to spend the vast majority of its time subdividing a large, empty, infeasible space just to find the two tiny islands where solutions exist.
The dReal baseline method might be more efficient at handling these kinds of disjoint feasible regions,
which would explain why our method is comparatively slower on this specific problem.

#subfigure(
  figure(
    image("../renders/full_comparison_Positive_Islands.svg", width: 110%),
    caption: [Positive islands runtime performance.],
  ),
  <fig:problem_runtime_positive_islands>,
  figure(
    image("../renders/full_comparison_Singularity_Edge.svg", width: 110%),
    caption: [Singularity edge runtime performance],
  ),
  <fig:problem_runtime_singularity_edge>,
  columns: (1fr, 1fr),
  caption: [(a) On certain problems, dReal defeats all other methods with marginal error.
    Our method is much slower than the dReal baseline (normalized runtime is around 6-7x slower the baseline time) while finding a similar bound.
    This is also the case for numeric optimizers. SHGO and Dual Annealing get the wrong answer.
    (b) On other problems, our method is able to match standard numerical optimization performance.
  ],
  label: <fig:problem_runtime_dreal_outlier>,
)



#figure(
  image("../renders/aggregate_normalized_runtime.svg", width: 110%),
  caption: [Plot of aggregate runtime performance across testing suite.],
) <fig:aggregate_runtime>

#figure(
  image("../renders/full_aggregate_normalized_runtime.svg", width: 110%),
  caption: [Plot of aggregate runtime performance across testing suite with standard numeric optimization methods.],
) <fig:aggregate_runtime_numeric>


// todo: make these subfigures

#figure(
  image("../renders/full_aggregate_bound_difference.svg", width: 110%),
  caption: [Plot of aggregate bound difference across testing suite.
    Note that the gradient-split and affine-bound heursitic combination has been ommitted as it tended to produce much larger (and incorrect) bounds on many problems.
    Dual annealing has been ommitted for a similar reason.],
) <fig:aggregate_bound_numeric>


//TODO: would be great to get visualizations on each of the particular problem plots

affine bounding heuristic can fail quickly
#figure(
  image("../renders/full_comparison_Split_Islands.svg", width: 60%),
  caption: [Plot of ],
) <fig:problem_runtime_split_islands>

Numerical methods sometimes outperform solver-aided ones
#subfigure(
  figure(image("../renders/full_comparison_Rational_Valley.svg"), caption: [Rational valley runtime performance.]),
  <fig:problem_runtime_rational_valley>,
  figure(
    image("../renders/full_comparison_Sanity_Poly.svg"),
    caption: [Sanity poly runtime performance],
  ),
  <fig:problem_runtime_sanity_poly>,
  columns: (1fr, 1fr),
  caption: [(a) On certain problems, dReal defeats all other methods with marginal error.
    (b) On other problems, our method is able to match standard numerical optimization performance.
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


