# TODO: extract algorithm errors to common constants

CONVERGENCE_TOLERANCE = 1e-10
MAX_STAGNANT_ITERS = 100

ERROR_INFEASIBLE = """
Found no feasible point for the current box under the given constraint.
"""

ERROR_NO_GLOBAL_MIN = """
Found no global minimum for the current objective over the given region.
"""
