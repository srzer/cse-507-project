
date around 11/20

Final approximate global lower bound B = 0.42072317943121723
Branch-and-bound time: 2.078697681427002 seconds

[baseline] approximate global minimum f ≈ 0.42098920842414606
Baseline dReal Minimize time: 2.0643882751464844 seconds


// after some organization & adding, it seems that i've managed to slow everything way down! uh oh.

objective function:
((3 - 2 * x0 - 2 * x1 - 2 * x2 + pow(x0, 2) + pow(x1, 2) + pow(x2, 2)) / (pow(x0, 2) + pow(x1, 2) + pow(x2, 2)))

initial feasible point:   0.525302882627382

Final approximate global lower bound B = 0.42072317943121723
Branch-and-bound time: 7.800956964492798 seconds

[baseline] approximate global minimum f ≈ 0.42094964010138736
Baseline dReal Minimize time: 9.014018535614014 seconds

// there seems to be a good deal of random variation in run times
Running with parameters:
  dim = 3
  min_box_size = 0.1
  delta = 0.001
  err = 0.0001
  init_box range: 1 to 10

=== Running GlobalMinBranchAndBound ===
Final approximate global lower bound B = 0.4207031351608186
Branch-and-bound time: 18.623109817504883 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
Improved Branch-and-bound time: 11.29397702217102 seconds

=== Running BaselineMin ===
[baseline] approximate global minimum f ≈ 0.42094964010138736
Baseline dReal Minimize time: 9.076762914657593 seconds

=== Running GlobalMinBranchAndBound ===
Final approximate global lower bound B = 0.4207031351608186
Branch-and-bound time: 17.57270121574402 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
Improved Branch-and-bound time: 11.118875980377197 seconds

=== Running BaselineMin ===
[baseline] approximate global minimum f ≈ 0.42094964010138736
Baseline dReal Minimize time: 8.472108125686646 seconds

=== Running GlobalMinBranchAndBound ===
Result: Right(value=0.4207031351608186)
Branch-and-bound time: 18.53873610496521 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
Result: Right(value=0.4207031351608186)
Improved Branch-and-bound time: 9.592769861221313 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 9.158195972442627 seconds
