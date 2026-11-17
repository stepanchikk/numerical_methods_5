"""Центральні похідні Стірлінґа. Частина Сандюка Даниїла."""
from .common import differences, uniform_step


def derivatives(x, y):
    if len(x) != 5 or len(y) != 5:
        raise ValueError("Ця формула потребує п'яти вузлів")
    h = uniform_step(x)
    d = differences(y)
    first = ((d[1][1]+d[1][2])/2
             - (d[3][0]+d[3][1])/12) / h
    second = (d[2][1] - d[4][0]/12) / h**2
    return first, second
