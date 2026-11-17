"""Перевірка точності на поліномах та незалежним інтегралом."""
import unittest
from src import newton, stirling, simple
from src.common import X, Y, A, B, integrand, exact_integral
from src.integration import quadrature, runge
from src.error_bounds import bounds


class MethodsTest(unittest.TestCase):
    def test_newton_exact_quartic(self):
        x = (-1.0, -0.5, 0.0, 0.5, 1.0)
        y = [v**4-2*v**3+3*v*v+v-7 for v in x]
        for backward in (False, True):
            for v in (-1, -0.3, 0.7, 1):
                first, second = newton.derivatives(x, y, v, backward)
                self.assertAlmostEqual(first, 4*v**3-6*v*v+6*v+1)
                self.assertAlmostEqual(second, 12*v*v-12*v+6)

    def test_stirling_and_newton(self):
        actual = stirling.derivatives(X, Y)
        expected = newton.derivatives(X, Y, X[2])
        for a, b in zip(actual, expected):
            self.assertAlmostEqual(a, b, places=12)

    def test_simple_quadratic(self):
        x = (0, 0.5, 1)
        self.assertEqual(simple.derivatives(x, [1, 1.75, 3], 1), (2.0, 2.0))

    def test_quadrature_polynomials(self):
        for method in ("midpoint", "trapezoid", "simpson"):
            self.assertAlmostEqual(quadrature(lambda x: 3*x+2, 0, 1, 6, method), 3.5)
        self.assertAlmostEqual(quadrature(lambda x: x**3, 0, 1, 6, "simpson"), 0.25)

    def test_runge_and_apriori(self):
        ref = exact_integral()
        limits = bounds(A, B)
        for method, eps, n in (("midpoint", .01, 5),
                                ("trapezoid", .001, 5), ("simpson", .0001, 6)):
            rows = runge(integrand, A, B, method, eps, n)
            self.assertLessEqual(rows[-1]["runge"], eps)
            self.assertLess(abs(rows[-1]["value"]-ref), eps)
            self.assertLessEqual(abs(quadrature(integrand, A, B, n, method)-ref), limits[method])

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            quadrature(integrand, A, B, 5, "simpson")
        with self.assertRaises(ValueError):
            newton.derivatives([0, 1, 3], [0, 1, 9], 1)


if __name__ == "__main__":
    unittest.main()
