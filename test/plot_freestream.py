import csv
import math
from pathlib import Path


def read_nodes(path: Path):
    xs, ys = [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
    return xs, ys


def read_elems(path: Path):
    tris = []
    r_l2 = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            tris.append((int(row["v0"]), int(row["v1"]), int(row["v2"])))
            r_l2.append(float(row["R_L2"]))
    return tris, r_l2


def read_solution(path: Path):
    rho, u, v, p = [], [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rho.append(float(row["rho"]))
            u.append(float(row["u"]))
            v.append(float(row["v"]))
            p.append(float(row["p"]))
    return rho, u, v, p


def main():
    nodes_path = Path("freestream_nodes.csv")
    elems_path = Path("freestream_elems.csv")
    sol_path = Path("freestream_solution_elems.csv")
    if not nodes_path.exists() or not elems_path.exists() or not sol_path.exists():
        raise SystemExit(
            "Cannot find freestream_nodes.csv / freestream_elems.csv / freestream_solution_elems.csv.\n"
            "Run ./test_freestream from your build dir first."
        )

    xs, ys = read_nodes(nodes_path)
    tris, r_l2 = read_elems(elems_path)
    rho, u, v, p = read_solution(sol_path)
    rmax = max(r_l2) if r_l2 else 0.0
    print(f"Loaded {len(tris)} triangles. max R_L2 = {rmax:.3e}")

    try:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except Exception as e:
        raise SystemExit(
            "matplotlib is required for plotting.\n"
            "Install: python3 -m pip install matplotlib\n"
            f"Import error: {e}"
        )

    tri = mtri.Triangulation(xs, ys, triangles=tris)
    logR = [math.log10(r + 1e-300) for r in r_l2]

    fig, axs = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    axs = axs.ravel()

    def panel(ax, facecolors, title, cbar_label):
        tpc = ax.tripcolor(tri, facecolors=facecolors, shading="flat", cmap="viridis")
        ax.set_aspect("equal", "box")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(tpc, ax=ax, label=cbar_label)

    panel(axs[0], logR, "log10(R_L2)", "log10(R_L2)")
    panel(axs[1], rho, "rho", "rho")
    panel(axs[2], u, "u", "u")
    panel(axs[3], v, "v", "v")
    panel(axs[4], p, "p", "p")

    axs[5].axis("off")
    fig.suptitle("Freestream on [0,1]x[0,1]: residual + solution fields")
    plt.show()


if __name__ == "__main__":
    main()

