"""Похідні Ньютона. Частина Загаєвського Олександра."""
from .common import differences, uniform_step


def derivatives(x, y, point, backward=False):
    """Перші дві похідні повного многочлена Ньютона."""
    if len(x) != len(y):
        raise ValueError("Кількість x і y має збігатися")
    h = uniform_step(x)
    table = differences(y)
    t = (point - (x[-1] if backward else x[0])) / h
    # Базис b та його похідні db, ddb за t.
    b, db, ddb = 1.0, 0.0, 0.0
    first, second = 0.0, 0.0
    for k in range(1, len(y)):
        factor = t + (k-1) if backward else t - (k-1)
        b, db, ddb = (b*factor/k, (db*factor+b)/k,
                       (ddb*factor+2*db)/k)
        delta = table[k][-1] if backward else table[k][0]
        first += delta * db
        second += delta * ddb
    return first/h, second/(h*h)
