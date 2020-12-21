import numpy as np


def f0(x, y):
    return -y - x * x


def f0_solution(x):
    return -x * x + 2 * x - 2 + 12 * (np.e ** (-x))


def f1(x, y):
    return y + 2 * x - 3


def f1_solution(x):
    return 1 + np.e ** x - 2 * x


def system1_y1(x, y):
    u = y[0]
    v = y[1]
    return (x * x + 0.1 * u * u) ** .5 + v


def system1_y2(x, y):
    u = y[0]
    v = y[1]
    return np.cos(2.1 * v) + u


def system2_y1(x, y):
    u = y[0]
    v = y[1]
    return np.cos(u + 1.1 * v) + 2.1


def system2_y2(x, y):
    u = y[0]
    v = y[1]
    return 1.1 / (x + 2.1 * (u * u)) + x + 1


def system3_y1(x, y):
    u = y[0]
    v = y[1]
    return -v


def system3_y2(x, y):
    u = y[0]
    v = y[1]
    return u


def system4_y1(x, y):
    u = y[0]
    v = y[1]
    return 3 * u - v


def system4_y2(x, y):
    u = y[0]
    v = y[1]
    return 4 * u - v


def system4_y1_solution(x):
    return (1 + x) * np.e ** x


def system4_y2_solution(x):
    return (1 + 2 * x) * np.e ** x


def system3_y1_solution(x):
    return np.cos(x)


def system3_y2_solution(x):
    return np.sin(x)


def p1(x):
    return -3 * x


def q1(x):
    return 2


def func1(x):
    return -1.5


def Bundle2_solution(x):
    return -2 + np.e ** x


def p2(x):
    return -1


def q2(x):
    return 0


def func2(x):
    return 0


def Bundle3_solution(x):
    return 1 - np.sin(x) - np.cos(x)


def p3(x):
    return 0


def q3(x):
    return 1


def func3(x):
    return 1

