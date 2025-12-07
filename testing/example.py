from dreal import And, Variable
from box import build_basic_box
from objective import Polynomial, Rational, Term, to_term
from .test_type import Problem

def make_example_problem_1(dim: int) -> Problem:
    """
    Creates the first example test problem.
    This problem was previously used as the main demonstration.
    """
    # --- Variables ---
    variables = [Variable(f"x{i}") for i in range(dim)]

    # --- Objective Function ---
    numerator_terms: list[Term] = []
    denominator_terms: list[Term] = []
    
    # Denominator: sum_i xi^2
    for i in range(dim):
        exps = [0] * dim
        exps[i] = 2
        denominator_terms.append(to_term(1.0, exps))

    # Numerator: sum_i (xi^2 - 2xi + 1)
    constant_term = 0.0
    for i in range(dim):
        exps2 = [0] * dim; exps2[i] = 2
        numerator_terms.append(to_term(1.0, exps2))
        exps1 = [0] * dim; exps1[i] = 1
        numerator_terms.append(to_term(-2.0, exps1))
        constant_term += 1.0
    numerator_terms.append(to_term(constant_term, [0] * dim))
    
    objective = Rational(Polynomial(numerator_terms), Polynomial(denominator_terms))

    # --- Constraints ---
    center1 = [3] * dim
    center2 = [4] * dim
    center3 = [3] + [4] * (dim - 1)
    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))

    constraints = And(
        sq_dist(variables, center1) <= 4,
        sq_dist(variables, center2) <= 4,
        sq_dist(variables, center3) <= 4,
    )

    # --- Initial Box ---
    initial_box = build_basic_box(1, 10, dim)

    return Problem(
        name="Example Problem 1",
        dim=dim,
        objective=objective,
        constraints=constraints,
        initial_box=initial_box,
        variables=variables
    )

def make_example_problem_2(dim: int) -> Problem:
    """
    Creates the second example test problem.
    """
    # --- Variables ---
    variables = [Variable(f"x{i}") for i in range(dim)]

    # --- Objective Function ---
    numerator_terms: list[Term] = []
    denominator_terms: list[Term] = []

    # Numerator
    for i in range(dim):
        exps = [0] * dim; exps[i] = 4
        numerator_terms.append(to_term(1.0, exps))
    for i in range(dim):
        exps = [0] * dim; exps[i] = 2
        numerator_terms.append(to_term(-3.0, exps))
    for i in range(dim):
        for j in range(i + 1, dim):
            exps = [0] * dim; exps[i] = 2; exps[j] = 2
            numerator_terms.append(to_term(0.5, exps))
    numerator_terms.append(to_term(2.0, [0] * dim))

    # Denominator
    denominator_terms.append(to_term(1.0, [0] * dim))
    for i in range(dim):
        exps = [0] * dim; exps[i] = 2
        denominator_terms.append(to_term(1.0, exps))

    objective = Rational(Polynomial(numerator_terms), Polynomial(denominator_terms))

    # --- Constraints ---
    center1 = [0.0] * dim
    center2 = [0.7] * dim
    center3 = [-0.7] * dim
    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))
    
    constraints = And(
        sq_dist(variables, center1) <= 16,
        sq_dist(variables, center2) <= 16,
        sq_dist(variables, center3) <= 16,
    )

    # --- Initial Box ---
    initial_box = build_basic_box(-10, 10, dim)

    return Problem(
        name="Example Problem 2",
        dim=dim,
        objective=objective,
        constraints=constraints,
        initial_box=initial_box,
        variables=variables
    )
