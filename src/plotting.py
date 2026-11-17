"""Графіки лабораторної. Спільна частина Шудровського Степана."""
from pathlib import Path
from .common import X, Y, A, B, differences, integrand, uniform_step
from .newton import derivatives


def interpolation(point):
    """Значення того ж многочлена Ньютона, що диференціюється."""
    t = (point-X[0])/uniform_step(X)
    table = differences(Y)
    value, basis = Y[0], 1.0
    for k in range(1, len(X)):
        basis *= (t-k+1)/k
        value += table[k][0]*basis
    return value


def parabola(point, x0, x1, x2):
    """Локальна парабола Лагранжа для однієї панелі Сімпсона."""
    return (integrand(x0)*(point-x1)*(point-x2)/((x0-x1)*(x0-x2))
            + integrand(x1)*(point-x0)*(point-x2)/((x1-x0)*(x1-x2))
            + integrand(x2)*(point-x0)*(point-x1)/((x2-x0)*(x2-x1)))


def grid(a, b, n=400):
    return [a+(b-a)*i/n for i in range(n+1)]


def create_plots(output_dir, show=False):
    """Зберегти два PNG; show=True також відкриває вікна графіків."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    fig, ax = plt.subplots(figsize=(8.2, 4.8), layout="constrained")
    xs = grid(X[0], X[-1])
    ax.plot(xs, [interpolation(x) for x in xs], color="#1e5a91",
            linewidth=2, label="Інтерполяційний многочлен P₄(x)")
    ax.scatter(X, Y, color="#b33235", s=45, zorder=4,
               label="Табличні значення")
    slope, _ = derivatives(X, Y, X[2])
    tangent_x = grid(X[2]-.4, X[2]+.4, 30)
    ax.plot(tangent_x, [Y[2]+slope*(x-X[2]) for x in tangent_x],
            color="#b67b12", linestyle="--", linewidth=1.8,
            label="Дотична в x = 1,531; P₄′ = −2,270")
    for i, (x, y) in enumerate(zip(X, Y)):
        ax.annotate(f"x{i}", (x, y), xytext=(7, 8),
                    textcoords="offset points", fontsize=10)
    ax.set(xlabel="x", ylabel="y", title="Таблична функція та її інтерполяція")
    ax.grid(alpha=.25)
    ax.legend(loc="lower left", fontsize=9)
    path = output_dir / "interpolation.png"
    fig.savefig(path, dpi=220)
    paths.append(path)

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 10.0),
                             layout="constrained", sharex=True)
    curve_x = grid(A, B)
    titles = ("Середні прямокутники, n = 5",
              "Трапеції, n = 5", "Сімпсон, n = 6")
    colors = ("#477eaf", "#ba782e", "#3b8b67")
    for ax, title, color in zip(axes, titles, colors):
        ax.set(title=title, ylabel="f(x)", ylim=(0, 1.5))
        ax.grid(alpha=.2)
        ax.set_axisbelow(True)
    h = (B-A)/5
    for i in range(5):
        left, right = A+i*h, A+(i+1)*h
        mid = (left+right)/2
        axes[0].add_patch(Rectangle((left, 0), h, integrand(mid),
                          facecolor=colors[0], alpha=.25,
                          edgecolor=colors[0], linewidth=1.5))
        axes[0].plot(mid, integrand(mid), "o", color=colors[0], ms=4)
        axes[1].add_patch(Polygon([(left, 0), (left, integrand(left)),
                                  (right, integrand(right)), (right, 0)],
                          facecolor=colors[1], alpha=.25,
                          edgecolor=colors[1], linewidth=1.5))
    h = (B-A)/6
    for i in range(0, 6, 2):
        x0, x1, x2 = A+i*h, A+(i+1)*h, A+(i+2)*h
        panel_x = grid(x0, x2, 60)
        panel_y = [parabola(x, x0, x1, x2) for x in panel_x]
        axes[2].fill_between(panel_x, 0, panel_y,
                             color=colors[2], alpha=.2+(i//2)*.07)
        axes[2].plot(panel_x, panel_y, color=colors[2], lw=2)
        axes[2].vlines([x0, x2], 0, [integrand(x0), integrand(x2)],
                       colors=colors[2], linewidth=1)
    for ax, color in zip(axes, colors):
        ax.plot(curve_x, [integrand(x) for x in curve_x],
                color="#222222", linestyle="--", linewidth=1.5,
                label="f(x) = (x + 2)/(x² − x + 1)")
        ax.plot([], [], color=color, linewidth=5, alpha=.5,
                label="Квадратурне наближення площі")
        ax.legend(loc="lower left", fontsize=9)
        ax.set_xlim(A, B)
    axes[-1].set_xlabel("x")
    path = output_dir / "integration_methods.png"
    fig.savefig(path, dpi=220)
    paths.append(path)
    if show:
        plt.show()
    plt.close("all")
    return paths
