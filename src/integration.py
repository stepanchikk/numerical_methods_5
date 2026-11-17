"""Квадратури та Рунге. Частина Тимчука Тараса."""
from math import fsum


def quadrature(f, a, b, n, method):
    if not isinstance(n, int) or n <= 0 or b <= a:
        raise ValueError("Потрібні n > 0 та a < b")
    h = (b-a)/n
    if method == "midpoint":
        return h*fsum(f(a+(i+0.5)*h) for i in range(n))
    if method == "trapezoid":
        return h*((f(a)+f(b))/2
                  + fsum(f(a+i*h) for i in range(1, n)))
    if method == "simpson":
        if n % 2:
            raise ValueError("Для Сімпсона n має бути парним")
        return h/3*(f(a)+f(b)
                    + 4*fsum(f(a+i*h) for i in range(1, n, 2))
                    + 2*fsum(f(a+i*h) for i in range(2, n, 2)))
    raise ValueError("Невідомий метод")


def runge(f, a, b, method, eps, n0, max_steps=20):
    if eps <= 0:
        raise ValueError("Точність має бути додатною")
    p = 4 if method == "simpson" else 2
    n = n0
    coarse = quadrature(f, a, b, n, method)
    rows = [{"n": n, "h": (b-a)/n, "value": coarse, "runge": None}]
    for _ in range(max_steps):
        n *= 2
        fine = quadrature(f, a, b, n, method)
        error = abs(fine-coarse)/(2**p-1)
        rows.append({"n": n, "h": (b-a)/n,
                     "value": fine, "runge": error})
        if error <= eps:
            return rows
        coarse = fine
    raise RuntimeError("Ліміт подвоєнь вичерпано")
