"""Спільний запуск. Тімлід — Шудровський Степан."""
import json
import argparse
from pathlib import Path
from src.common import X, Y, A, B, differences, exact_integral, integrand
from src import newton, stirling, simple
from src.integration import runge, quadrature
from src.error_bounds import bounds

METHODS = {"midpoint": "Середні прямокутники",
           "trapezoid": "Трапеції", "simpson": "Сімпсон"}


def calculate():
    derivatives = []
    for point, backward in ((0.231, False), (2.2, True)):
        d1, d2 = newton.derivatives(X, Y, point, backward)
        derivatives.append({"method": "Ньютон", "x": point,
                            "first": d1, "second": d2})
    d1, d2 = stirling.derivatives(X, Y)
    derivatives.append({"method": "Стірлінґ", "x": X[2],
                        "first": d1, "second": d2})
    for i in (2, 3):
        d1, d2 = simple.derivatives(X, Y, i)
        derivatives.append({"method": "Три вузли", "x": X[i],
                            "first": d1, "second": d2})
    reference = exact_integral()
    apriori = bounds(A, B)
    integrals = {}
    for method, eps, n0 in (("midpoint", 1e-2, 5),
                            ("trapezoid", 1e-3, 5),
                            ("simpson", 1e-4, 6)):
        rows = runge(integrand, A, B, method, eps, n0)
        fixed = quadrature(integrand, A, B, n0, method)
        integrals[method] = {"eps": eps, "rows": rows,
                            "actual_error": abs(rows[-1]["value"]-reference),
                            "fixed_n": n0, "fixed_value": fixed,
                            "fixed_error": abs(fixed-reference)}
    return {"variant": 7, "x": X, "y": Y,
            "differences": differences(Y), "derivatives": derivatives,
            "exact_integral": reference, "apriori": apriori,
            "integrals": integrals}


def format_results(data):
    lines = ["Лабораторна робота №5. Варіант 7", "", "Таблиця різниць"]
    for k, row in enumerate(data["differences"]):
        lines.append(f"Δ^{k}: " + "  ".join(f"{v:.6f}" for v in row))
    lines.append("\nПохідні: метод | x | f′ | f″")
    for row in data["derivatives"]:
        lines.append(f"{row['method']} | {row['x']:.3f} | "
                     f"{row['first']:.9f} | {row['second']:.9f}")
    lines.append(f"\nКонтрольний інтеграл: {data['exact_integral']:.12f}")
    for method, part in data["integrals"].items():
        lines.append(f"\n{METHODS[method]}, ε={part['eps']}")
        for row in part["rows"]:
            error = "—" if row["runge"] is None else f"{row['runge']:.10g}"
            lines.append(f"n={row['n']}, h={row['h']:.9f}, "
                         f"I={row['value']:.12f}, Рунге={error}")
        lines.append(f"Фактична похибка: {part['actual_error']:.10g}")
        lines.append(f"Апріорна межа при n={part['fixed_n']}: "
                     f"{data['apriori'][method]:.10g}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="ЛР5, варіант 7")
    parser.add_argument("--show-plots", action="store_true",
                        help="Показати графіки у вікнах після збереження")
    args = parser.parse_args()
    data = calculate()
    text = format_results(data)
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    (output / "results.json").write_text(
        json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "results.txt").write_text(text + "\n", encoding="utf-8")
    print(text)
    from src.plotting import create_plots
    paths = create_plots(output / "plots", show=args.show_plots)
    print("\nГрафіки збережено:")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
