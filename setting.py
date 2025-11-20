from dataclasses import dataclass
from dreal import *
@dataclass
class Box3D:
    xl: float; xu: float
    yl: float; yu: float
    zl: float; zu: float

def f_value(x, y, z):
    return ((x-1)**2 + (y-1)**2 + (z-1)**2) / (x + y + z)

def f_constraint(x, y, z):
    return And((x-3)**2 + (y-3)**2 + (z-3)**2 <= 2 * 2, (x-4)**2 + (y-4)**2 + (z-4)**2 <= 2 * 2, (x-3)**2 + (y-4)**2 + (z-4)**2 <= 2 * 2)
    # return x * x + y * y + z * z <= 2 * 2

def get_init_box() -> Box3D:
    return Box3D(1, 10, 1, 10, 1, 10), 0.1

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

def get_init_box() -> Box3D:
    # huge outer box, feasible region is extremely sparse inside
    return Box3D(-50, 50, -50, 50, -50, 50), 2.0