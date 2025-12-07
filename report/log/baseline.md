
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

=== Running GlobalMinBranchAndBound ===
  Iteration 100, queue size: 6, current lower bound: 0.44155163744897263
  Iteration 200, queue size: 2, current lower bound: 0.43136736935460956
  Iteration 300, queue size: 4, current lower bound: 0.4309169917600907
  Iteration 400, queue size: 4, current lower bound: 0.4309169917600907
  Iteration 500, queue size: 8, current lower bound: 0.4216352618445757
algorithm completed after 579 iterations
Result: Right(value=0.4207031351608186)
Branch-and-bound time: 18.555098295211792 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
Result: Right(value=0.4207031351608186)
Improved Branch-and-bound time: 9.594598293304443 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 9.063214302062988 seconds

=== Running GlobalMinBranchAndBound ===
  iteration 100, queue size: 6, current lower bound: 0.44155163744897263
  iteration 200, queue size: 2, current lower bound: 0.43136736935460956
  iteration 300, queue size: 4, current lower bound: 0.4309169917600907
  iteration 400, queue size: 4, current lower bound: 0.4309169917600907
  iteration 500, queue size: 8, current lower bound: 0.4216352618445757
algorithm completed after 579 iterations
Result: Right(value=0.4207031351608186)
Branch-and-bound time: 17.61022686958313 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
Result: Right(value=0.4207031351608186)
Improved Branch-and-bound time: 9.167948007583618 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 8.539322137832642 seconds

// changing up the splitting and bounding techniques messes things up it seems
running with parameters:
  dim          = 3
  init_box     = BoxN(min=Point(coords=(1, 1, 1)), max=Point(coords=(10, 10, 10)))
  obj fn       = Rational(num=Polynomial(terms=[(1.0, (2, 0, 0)), (-2.0, (1, 0, 0)), (1.0, (0, 2, 0)), (-2.0, (0, 1, 0)), (1.0, (0, 0, 2)), (-2.0, (0, 0, 1)), (3.0, (0, 0, 0))]), den=Polynomial(terms=[(1.0, (2, 0, 0)), (1.0, (0, 2, 0)), (1.0, (0, 0, 2))]))
  contraint    = (((pow((-4 + x0), 2) + pow((-4 + x1), 2) + pow((-4 + x2), 2)) <= 4) and ((pow((-4 + x1), 2) + pow((-4 + x2), 2) + pow((-3 + x0), 2)) <= 4) and ((pow((-3 + x0), 2) + pow((-3 + x1), 2) + pow((-3 + x2), 2)) <= 4))
  min_box_size = 0.1
  delta        = 0.001
  err          = 0.0001
  splitter     = <box.split.SplitGradient object at 0x7f6d39dc8910>
  bounder      = <objective.bound.affine.AffineBounds object at 0x7f6d39dcb730>

=== Running GlobalMinBranchAndBound ===                                                                                                                                                                     iteration 100, queue size: 56, current lower bound: 0.525302882627382, stagnant: 100
  iteration 200, queue size: 156, current lower bound: 0.525302882627382, stagnant: 200                                                                                                                     iteration 300, queue size: 256, current lower bound: 0.525302882627382, stagnant: 300
  iteration 400, queue size: 356, current lower bound: 0.525302882627382, stagnant: 400
  iteration 500, queue size: 456, current lower bound: 0.525302882627382, stagnant: 500
  iteration 600, queue size: 556, current lower bound: 0.525302882627382, stagnant: 600
  iteration 700, queue size: 656, current lower bound: 0.525302882627382, stagnant: 700
  iteration 800, queue size: 756, current lower bound: 0.525302882627382, stagnant: 800
  iteration 900, queue size: 856, current lower bound: 0.525302882627382, stagnant: 900
  iteration 1000, queue size: 956, current lower bound: 0.525302882627382, stagnant: 1000
  converged: no improvement for 1000 iterations
  final lower bound: 0.525302882627382
algorithm completed after 1001 iterations
Result: Right(value=0.525302882627382)
Branch-and-bound time: 10.757195949554443 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
  iteration 100, queue size: 56, current lower bound: 0.525302882627382, stagnant: 100                                                                                                                      iteration 200, queue size: 156, current lower bound: 0.525302882627382, stagnant: 200
  iteration 300, queue size: 256, current lower bound: 0.525302882627382, stagnant: 300                                                                                                                     iteration 400, queue size: 356, current lower bound: 0.525302882627382, stagnant: 400
  iteration 500, queue size: 456, current lower bound: 0.525302882627382, stagnant: 500
  iteration 600, queue size: 556, current lower bound: 0.525302882627382, stagnant: 600
  iteration 700, queue size: 656, current lower bound: 0.525302882627382, stagnant: 700
  iteration 800, queue size: 756, current lower bound: 0.525302882627382, stagnant: 800
  iteration 900, queue size: 856, current lower bound: 0.525302882627382, stagnant: 900
  iteration 1000, queue size: 956, current lower bound: 0.525302882627382, stagnant: 1000
  Converged: no improvement for 1000 iterations
  Final lower bound: 0.525302882627382
Result: Right(value=0.525302882627382)
Improved Branch-and-bound time: 10.829509735107422 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 9.055834770202637 seconds

// longest side seems to be good
running with parameters:
  dim          = 3
  init_box     = BoxN(min=Point(coords=(1, 1, 1)), max=Point(coords=(10, 10, 10)))
  obj fn       = Rational(num=Polynomial(terms=[(1.0, (2, 0, 0)), (-2.0, (1, 0, 0)), (1.0, (0, 2, 0)), (-2.0, (0, 1, 0)), (1.0, (0, 0, 2)), (-2.0, (0, 0, 1)), (3.0, (0, 0, 0))]), den=Polynomial(terms=[(1.0, (2, 0, 0)), (1.0, (0, 2, 0)), (1.0, (0, 0, 2))]))
  contraint    = (((pow((-4 + x0), 2) + pow((-4 + x1), 2) + pow((-4 + x2), 2)) <= 4) and ((pow((-4 + x1), 2) + pow((-4 + x2), 2) + pow((-3 + x0), 2)) <= 4) and ((pow((-3 + x0), 2) + pow((-3 + x1), 2) + pow((-3 + x2), 2)) <= 4))
  min_box_size = 0.1
  delta        = 0.001
  err          = 0.0001
  splitter     = <box.split.SplitLongestSide object at 0x7f9d6f2cbd30>
  bounder      = <objective.bound.affine.AffineBounds object at 0x7f9d6f2cac50>

=== Running GlobalMinBranchAndBound ===
  iteration 100, queue size: 6, current lower bound: 0.44155163744897263, stagnant: 8
  iteration 200, queue size: 2, current lower bound: 0.43136736935460956, stagnant: 42
  iteration 300, queue size: 4, current lower bound: 0.4309169917600907, stagnant: 24
  converged: no improvement for 100 iterations
  final lower bound: 0.4309169917600907
algorithm completed after 377 iterations
Result: Right(value=0.4309169917600907)
Branch-and-bound time: 12.886556625366211 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
  iteration 100, queue size: 6, current lower bound: 0.4480440308326043, stagnant: 16
  iteration 200, queue size: 2, current lower bound: 0.43136736935460956, stagnant: 42
  iteration 300, queue size: 4, current lower bound: 0.4309169917600907, stagnant: 24
  Converged: no improvement for 100 iterations
  Final lower bound: 0.4309169917600907
Result: Right(value=0.4309169917600907)
Improved Branch-and-bound time: 3.3941075801849365 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 9.033820629119873 seconds

// bernstein is fast but bad??
running with parameters:
  dim          = 3
  init_box     = BoxN(min=Point(coords=(1, 1, 1)), max=Point(coords=(10, 10, 10)))
  obj fn       = Rational(num=Polynomial(terms=[(1.0, (2, 0, 0)), (-2.0, (1, 0, 0)), (1.0, (0, 2, 0)), (-2.0, (0, 1, 0)), (1.0, (0, 0, 2)), (-2.0, (0, 0, 1)), (3.0, (0, 0, 0))]), den=Polynomial(terms=[(1.0, (2, 0, 0)), (1.0, (0, 2, 0)), (1.0, (0, 0, 2))]))
  contraint    = (((pow((-4 + x0), 2) + pow((-4 + x1), 2) + pow((-4 + x2), 2)) <= 4) and ((pow((-4 + x1), 2) + pow((-4 + x2), 2) + pow((-3 + x0), 2)) <= 4) and ((pow((-3 + x0), 2) + pow((-3 + x1), 2) + pow((-3 + x2), 2)) <= 4))
  min_box_size = 0.1
  delta        = 0.001
  err          = 0.0001
  splitter     = <box.split.SplitGradient object at 0x7f19f04cac50>
  bounder      = <objective.bound.bernstein.BernsteinBounds object at 0x7f19f04cbd30>

=== Running GlobalMinBranchAndBound ===
  iteration 100, queue size: 56, current lower bound: 0.525302882627382, stagnant: 100
  converged: no improvement for 100 iterations
  final lower bound: 0.525302882627382
algorithm completed after 101 iterations
Result: Right(value=0.525302882627382)
Branch-and-bound time: 1.1723880767822266 seconds

=== Running ImprovedGlobalMinBranchAndBound ===
  iteration 100, queue size: 56, current lower bound: 0.525302882627382, stagnant: 100
  Converged: no improvement for 100 iterations
  Final lower bound: 0.525302882627382
Result: Right(value=0.525302882627382)
Improved Branch-and-bound time: 1.2103965282440186 seconds

=== Running BaselineMin ===
Result: Right(value=0.42094964010138736)
Baseline dReal Minimize time: 8.513031005859375 seconds

// affine + gradient seems to not work
// but affine by itself seems to do pretty good!
