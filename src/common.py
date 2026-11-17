"""Спільні дані. Частина Шудровського Степана."""
from math import atan, log, sqrt, isclose

X = (0.231, 0.881, 1.531, 2.181, 2.831)
Y = (-2.748, -3.225, -3.898, -5.908, -6.506)
A, B = 2.0, 2.5


def uniform_step(x):
    if len(x) < 2:
        raise ValueError("Потрібно щонайменше два вузли")
    h = x[1] - x[0]
    if h <= 0 or any(not isclose(x[i+1]-x[i], h,
                                 rel_tol=1e-10, abs_tol=1e-12)
                     for i in range(len(x)-1)):
        raise ValueError("Потрібні зростаючі рівновіддалені вузли")
    return h


def differences(y):
    table = [list(y)]
    while len(table[-1]) > 1:
        row = table[-1]
        table.append([row[i+1]-row[i] for i in range(len(row)-1)])
    return table


def integrand(x):
    return (x + 2) / (x*x - x + 1)


def primitive(x):
    return 0.5*log(x*x-x+1) + 5/sqrt(3)*atan((2*x-1)/sqrt(3))


def exact_integral(a=A, b=B):
    return primitive(b) - primitive(a)
