"""Апріорні оцінки. Частина Шудровського Степана."""
from fractions import Fraction
from math import ceil


def derivative_bounds():
    # f'' додатна та спадна на [2; 2.5]: M2=f''(2)=10/9.
    m2 = 10/9
    # u=x-2. N4(u) та D(u) додатні й зростають.
    # На [w; v] маємо f⁽⁴⁾ <= N4(v)/D(w)^5.
    # Fraction виключає похибки округлення при пошуку межі.
    upper = Fraction(0)
    for i in range(500):
        w, v = Fraction(i, 1000), Fraction(i+1, 1000)
        numerator = 216+2160*v+3600*v**2+2160*v**3+480*v**4+24*v**5
        denominator = (3+3*w+w*w)**5
        upper = max(upper, numerator/denominator)
    return m2, ceil(upper*100)/100


def bounds(a, b):
    if a != 2 or b != 2.5:
        raise ValueError("Межі похідних доведені лише для [2; 2.5]")
    m2, m4 = derivative_bounds()
    h5, h6 = (b-a)/5, (b-a)/6
    return {"M2": m2, "M4": m4,
            "midpoint": (b-a)*h5**2*m2/24,
            "trapezoid": (b-a)*h5**2*m2/12,
            "simpson": (b-a)*h6**4*m4/180}
