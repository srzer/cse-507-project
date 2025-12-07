import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from test_suite import (
    run_test_suite,
    test_sanity_poly,
    test_sanity_rational,
    test_rational_bowl,
    test_himmelblau_ratio,
    test_split_islands,
    test_singularity_edge,
    test_pole_avoidance,
    test_rational_valley,
    test_positive_islands,
    test_sparse_intersection,
    test_main_example
)

if __name__ == "__main__":
    problems = [
        test_sanity_poly(),
        test_sanity_rational(),
        test_rational_bowl(),
        # test_himmelblau_ratio(),
        test_split_islands(),
        test_singularity_edge(),
        test_pole_avoidance(),
        test_rational_valley(),
        test_positive_islands(),
        test_sparse_intersection(),
        test_main_example()
    ]

    run_test_suite(problems)
