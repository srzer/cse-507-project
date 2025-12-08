# UW CSE507 Course Project

Optimizing positive rational functions $\R^n -> \R$ with polynomial inequality constraints using interval arithmetic, Bernstein polynomial form, and branch-and-bound methods.

## Current implementation of n-dimensional branch-and-prune algorithm

- The input rational function is represented as two lists of terms, each term is [coefficient, [exponent_vector]].
For example, $f(x,y) = (2x^2y + 3y^2) / (x^2 + y)$, then numerator_terms = [ [2, [2,1]], [3, [0,2]] ], and denominator_terms = [ [1, [2,0]], [1, [0,1]] ].

- We start with an initial box that is assumed to cover all feasible areas and contain no poles.
We also define a minimal box size, which is, by default, the as same as our error bound $\delta$ as used in dReal (e.g. 1e-3).

- For each box, the pipeline is:

    - Solve (Bernstein method) a rough lower bound on this box. If it's larger than current lower bound, skip.

    - Solve (dReal) a solution in the box satisfying the constraints and improving the lower bound at least by $\epsilon$.
        - If unsat, then skip;
        - If sat, then update the lower bound using the solution.
    - If the box is already smaller than the minimal box size, skip.

    - Check (dReal) whether the box is fully feasible.
        - If so, we invoke a sub-program (no dReal involved) to solve the minimum of this feasible box;
          - Solve (Bernstein method) a rough lower bound on this box. If it's larger than current lower bound, skip.
          - If the box is already smaller than the minimal box size, skip.
          - Otherwise split it into two small boxes and iterate.
        - Otherwise, split it into two small boxes.

### Algorithm Heuristics

Box feasibility sampling heuristics
- completely sound
- grid (with some density) (WIP)
- random (with some seed?) (WIP)

Box slicing heuristics 
- bisect on random dimension
- bisect on longest side dimension
- dimension of greatest change (derivative)

### Logging

The main script logs the runtime and lower bound of the algorithms on the test suite.

## Known Limitations

### dReal Solver Hangs

Some complex problem configurations, particularly those involving intricate constraints or objective functions with singularities, can cause the underlying `dreal` SMT solver to enter a very long-running computation, effectively "hanging" the algorithm.

This is a fundamental limitation of using an external solver without a direct timeout API. Due to the way Python interacts with the `dreal` C++ library (specifically, an issue with serializing `dreal` objects for inter-process communication), a robust, low-level timeout on individual `dreal` calls is not feasible in the current architecture.

To manage this, a skip list has been implemented in `main.py`:

```python
SKIP_LIST = [
    # ... configurations known to hang are listed here ...
]
```

Any test run that matches a configuration in this list will be skipped automatically, allowing the test suite to complete reliably. If you encounter a new, consistently hanging test, please add it to this list.


## To Do
[x] Establish a toolkit of smart (e.g. function-behavior-aware) splitting heuristics.

[x] Replace expensive dReal solver calls with the Bernstein or Affine methods.

[ ] Write more tests to demonstrate determine our method's advantage compared
to dReal.

[ ] Write variable spread heuristics (heuristic_splitting 240-283, 285-319) for smarter box splitting


## Requirements

- Python >=3.9
- dReal >4.20

### Setup dReal

(a) For local dReal installation, please refer to https://github.com/dreal/dreal4.

(b) To use the provided dReal container,
build the image with
```bash
docker compose build dreal
```
This will take a minute.

Run scripts with the containerized binary via
```bash
docker compose run --rm dreal main.py
```

Enter into the container shell with
```bash
docker compose run --rm --entrypoint bash dreal
```

This container includes:
- Ubuntu 22.04 base
- dReal (v4.21.06.2) binaries
- Python 3 with dreal and other packages in `requirements.txt`

#### dReal Python Package

The [dReal package](https://pypi.org/project/dreal/) may be installed with
```bash
pip install dreal
```
Note that the package does not provide some bindings e.g. `dreal.Formula` by default.
These can be provided by installing the required `.pyi` header[^1] to your environment, e.g. in `venv/lib/python3.*/site-packages/dreal/`.
```bash
pip install mypy
stubgen -p dreal -o stubs/
```

However, generating the header requires the dReal binary and will not work with the standalone `dreal` package from `pip`.
For this reason, we have included a copy of the (dReal v4.21.06.2) headers in [dreal-container/_dreal_py.pyi](./dreal-container/_dreal_py.pyi).


[^1]: The `__init__.py` header is already provided.
