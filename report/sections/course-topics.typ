Emergent discovery of poor solver performances guided our optimization choices.
While we learned some of the essential algorithms and techniques that have built up SMT tools, in practice, SMT tools can still feel like black boxes.
As mentioned in class, solvers are tightly coupled tools that exchange modularity for performance. A lack of modularity impairs understanding.
#footnote[We experienced this on the project level, with untyped code.]
Without a good understanding of dReal itself, we were optimizing slow SMT run times by reducing complexity of solver input, cutting down our search space using branch-and-bound techniques.
With a more explicit understanding of the dReal solver backend, we can both guide the construction of our heuristic simplification techniques, and,
even better, directly focus on the solver’s shortcomings instead of pre-processing our input.
That being said, we did encounter DPLL(ICP) as a case of the DPLL algorithm specialized to interval constraint propagation of SMT using interval arithmetic.
Certainly, it is hard to fit everything in during a quarter.
We can either focus more theoretically on how a solver backend works versus more practically on how to build solver aided tools.
The project provides an experiential way of learning the latter skill, however, as discussed above,
we didn’t quite have the knowledge or toolkit to unpack the SMT solver and build our optimization with respect to its internal implementation.
