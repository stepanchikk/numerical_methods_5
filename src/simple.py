"""Прості формули на трьох вузлах. Частина Сандюка Даниїла."""
from .common import uniform_step


def derivatives(x, y, index):
    if len(x) != len(y) or not 0 < index < len(y)-1:
        raise ValueError("Потрібен внутрішній вузол таблиці")
    h = uniform_step(x)
    first = (y[index+1] - y[index-1]) / (2*h)
    second = (y[index+1] - 2*y[index] + y[index-1]) / h**2
    return first, second
