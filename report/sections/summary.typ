#import "@preview/elsearticle:1.1.0": *

Below is Fig. @fig:aggregate_runtime.


#figure(
  image("../renders/aggregate_normalized_runtime.svg", width: 90%),
  caption: [Plot of aggregate runtime performance across testing suite.],
) <fig:aggregate_runtime>

#figure(
  image("../renders/full_aggregate_normalized_runtime.svg", width: 90%),
  caption: [Plot of aggregate runtime performance across testing suite with standard numeric optimization methods.],
) <fig:aggregate_runtime_numeric>


// todo: make these subfigures

#figure(
  image("../renders/full_aggregate_bound_difference.svg", width: 90%),
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

sometimes our methods worked well (because of specialization to rational functions?)
#figure(
  image("../renders/full_comparison_Pole_Avoidance.svg", width: 60%),
  caption: [Plot of ],
) <fig:problem_runtime_pole_avoidance>

dReal is sometimes an outlier
#subfigure(
  figure(image("../renders/full_comparison_Positive_Islands.svg"), caption: [Positive islands runtime performance.]),
  <fig:problem_runtime_positive_islands>,
  figure(
    image("../renders/full_comparison_Singularity_Edge.svg"),
    caption: [Singularity edge runtime performance],
  ),
  <fig:problem_runtime_singularity_edge>,
  columns: (1fr, 1fr),
  caption: [(a) On certain problems, dReal defeats all other methods with marginal error.
    (b) On other problems, our method is able to match standard numerical optimization performance.
  ],
  label: <fig:problem_runtime_dreal_outlier>,
)

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
== Picking Test Suite
