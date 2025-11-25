from dataclasses import dataclass
from dreal import *
from dataclasses import dataclass
from dreal import *

@dataclass
class BoxND:
    """Axis-aligned box in n dimensions."""
    lows: list
    highs: list

    def __post_init__(self):
        assert len(self.lows) == len(self.highs), "Dimension mismatch in BoxND."

    @property
    def dim(self):
        return len(self.lows)


def f_value(*xs):
    """
    Example objective function generalized to n dimensions:

        f(x) = sum_i (xi - 1)^2 / sum_i xi

    You can change this to your actual n-D function if needed.
    """
    num = sum((x - 1) ** 2 for x in xs)
    denom = sum(x**2 for x in xs)
    return num / denom


def f_constraint(*xs):
    """
    Example constraint generalized to n dimensions:

    The original 3D version had three ball constraints.
    Here we mimic a similar structure by placing several
    2-radius balls at fixed positions in n-D space.

    You can customize this if your real constraint is different.
    """

    # Build three centers similar to the 3D case:
    # (3,3,3,...), (4,4,4,...), (3,4,4,...)
    n = len(xs)

    center1 = [3] * n
    center2 = [4] * n
    center3 = [3] + [4] * (n - 1)

    def sq_dist(xs, center):
        return sum((x - c) ** 2 for x, c in zip(xs, center))

    return And(
        sq_dist(xs, center1) <= 2 * 2,
        sq_dist(xs, center2) <= 2 * 2,
        sq_dist(xs, center3) <= 2 * 2,
    )


def get_init_box(n: int):
    """
    Build an initial axis-aligned n-D box:

        [1,10] × [1,10] × ... × [1,10]

    Same as your original 3D version.
    """
    lows = [1.0] * n
    highs = [10.0] * n
    min_box_size = 0.5

    return BoxND(lows, highs), min_box_size

# @dataclass
# class Box3D:
#     xl: float; xu: float
#     yl: float; yu: float
#     zl: float; zu: float

# def f_value(x, y, z):
#     return ((x-1)**2 + (y-1)**2 + (z-1)**2) / (x + y + z)

# def f_constraint(x, y, z):
#     return And((x-3)**2 + (y-3)**2 + (z-3)**2 <= 2 * 2, (x-4)**2 + (y-4)**2 + (z-4)**2 <= 2 * 2, (x-3)**2 + (y-4)**2 + (z-4)**2 <= 2 * 2)
#     # return x * x + y * y + z * z <= 2 * 2

# def get_init_box() -> Box3D:
#     return Box3D(1, 10, 1, 10, 1, 10), 0.1

# def f_value(x, y, z):
#     # main well near (7, 3, 1)
#     num = (x - 8)**2 + (y - 3)**2 + (z - 1)**2 + 0.1 * ((x + 3)**2 + (y - 2)**2)
#     den = 1 + 0.02 * (x * x + y * y + z * z)
#     return num / den

# def f_constraint(x, y, z):
#     # many tiny balls (radius R) scattered in a huge box
#     R = 0.1

#     def ball(cx, cy, cz):
#         return (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= R * R

#     # 5 small balls in far apart locations
#     c1 = ball(8,  3,  1)
#     c2 = ball(20, 22, 20)
#     c3 = ball(-10, 26, -15)
#     c4 = ball(32, -5, -7)
#     c5 = ball(-26, -30, 10)

#     return Or(c1, Or(c2, Or(c3, Or(c4, c5))))

# def get_init_box() -> Box3D:
    # huge outer box, feasible region is extremely sparse inside
    # return Box3D(-50, 50, -50, 50, -50, 50), 2.0