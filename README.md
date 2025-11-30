# UW CSE507 Course Project

Optimizing positive rational functions $\R^n -> \R$ with polynomial inequality constraints using interval arithmetic, Bernstein polynomial form, and branch-and-bound methods.

## Current implementation of n-dimensional branch-and-prune algorithm using Bernstein method

- The input rational function is represented as two lists of terms, each term is [coefficient, [exponent_vector]].
For example, $f(x,y) = (2x^2y + 3y^2) / (x^2 + y)$, then numerator_terms = [ [2, [2,1]], [3, [0,2]] ], and denominator_terms = [ [1, [2,0]], [1, [0,1]] ].

- We start with an initial box that is assumed to cover all feasible areas and contain no poles.
We also define a minimal box size, which is, by default, the as same as our error bound $\delta$ as used in dReal (e.g. 1e-3).

- For each box, the pipeline is:

    - Solve (Bernstein method) a rough lower bound on this box. If it's larger than current lower bound, skip.

    - Solve (dReal) a solution in the box satisfying the constraints and improving the lower bound at least by $\epsilon$.
    If unsat, then skip; If sat, then update the lower bound.
    - If the box is already smaller than the minimal box size, skip.

    - Check (dReal) if the box is fully feasible.
    If so, we solve (dReal) the accurate lower bound with only box constraint; otherwise, we split the box (currently splitting the longest edge).

### Algorithm Heuristics

Box feasibility sampling heuristics
- completely sound
- grid (with some density)
- random (with some seed?)

Box slicing heuristics 
- 1?
- 2?
- dimension of greatest change (derivative)

## To Do
- Establish a toolkit of smart (e.g. function-behavior-aware) splitting heuristics.
- Replace expensive dReal solver calls through the Bernstein method.
- Write more tests to demonstrate our method's advantage.

## Requirements

- Python >=3.9
- dReal >4.20

### Setup dReal

(a) For local dReal installation, please refer to https://github.com/dreal/dreal4.

(b) To use the provided dReal container,
set the `ROOT_DIR` environment variable using `make`: [^1]
```bash
make .env
```

build the image with
```bash
docker compose build dreal
```

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
- Python 3 with dreal package

#### dReal Python Package

The [dReal package](https://pypi.org/project/dreal/) may be installed with
```bash
pip install dreal
```
Note that the package does not provide some bindings e.g. `dreal.Formula` by default.
These can be provided by installing the required `.pyi` header[^2] to your environment, e.g. in `venv/lib/python3.*/site-packages/dreal/`.
```bash
pip install mypy
stubgen -p dreal -o stubs/
```

However, generating the header requires the dReal binary and will not work with the standalone `dreal` package from `pip`.
For this reason, we have included a copy of the (dReal v4.21.06.2) headers in [dreal-container/_dreal_py.pyi](./dreal-container/_dreal_py.pyi).



#### Python dReal package

The [dReal package](https://pypi.org/project/dreal/) may be installed with
```bash
pip install dreal
```
Note that the package does not provide some bindings e.g. `dreal.Formula` by default.
We can provide these to some 
```bash
pip install mypy
stubgen -p dreal -o stubs/
```



```bash
pip install mypy
stubgen -p dreal -o stubs/
```


[^1]: Alternatively, without `make`, `cd` into the project root directory and run `echo ROOT_DIR=$(pwd) > .env`.
[^2]: The `__init__.py` header is already provided.
