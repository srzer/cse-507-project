Jaedon:
Researched and tested different bounding algorithms across various sample problems.
Explored Bernstein polynomials and full and partial affine arithmetic for pruning effectiveness, and considered different heuristics for box splitting:
splitting along the objective function gradient, splitting along the gradient of the constraints, splitting along the widest dimension, splitting along an aggregate of these heuristics.

Ruizhe:
Implemented the basic branch-and-bound algorithmic framework based on dReal, and conducted comparisons against dReal’s Minimize procedure.
Integrated the Bernstein method, the affine method and the gradient-guided splitting heuristic into the branch-and-bound framework.
Initially designed simple test cases to demonstrate the advantages of our proposed methods.

Evan:
Tested our implementations on various parameters and optimization functions to find strengths and weaknesses of our code.
Produced visualizations in Blender (supported by scripts) by extracting boxes and current best points to provide better intuition to support the results from the tests.

lipson:
Took meeting minutes for logistics: ensuring we had actionable tasks, planning the next meeting date, summarizing our current state and what we need to implement next.
Unified python scripts with types and interfaces to be on the same page and communicate about our code.
In the future, it could have been better to start with everyone on the same page about code quality practices, e.g. strong typing.
Kept project documentation up to date: how to use our tool, install dReal container.
