#!/usr/bin/env python3
"""
Generate a q=2 curved .gri mesh from a linear (q=1) triangular .gri mesh.

Workflow:
1) Read linear mesh (.gri), blade coordinate files (upper/lower).
2) Curve the blade-wall-adjacent elements (default: BGroup2 and BGroup6).
3) Write a mixed-order .gri:
   - curved elements:  p=2 TriLagrange (6 nodes)
   - remaining elems:  p=1 TriLagrange (3 nodes)
4) Save before/after visualization.

This script does not modify solver code.
"""

from __future__ import annotations

import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar


@dataclass
class BoundaryGroup:
    nface: int
    nnode: int
    title: str
    faces: List[List[int]]  # 1-based node ids from file


@dataclass
class ElementGroup:
    nelem: int
    p: int
    basis: str
    elems: List[List[int]]  # 1-based node ids from file


@dataclass
class PeriodicGroup:
    npair: int
    kind: str
    pairs: List[Tuple[int, int]]  # 1-based node ids from file


@dataclass
class GriMesh:
    nnode: int
    nelem: int
    dim: int
    vertices: List[Tuple[float, float]]
    boundaries: List[BoundaryGroup]
    element_groups: List[ElementGroup]
    periodic_groups: List[PeriodicGroup]


def _import_pyplot():
    # Avoid non-writable default cache warnings on some environments.
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = Path(tempfile.gettempdir()) / "matplotlib-cache"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    return plt


def parse_gri(path: Path) -> GriMesh:
    lines = path.read_text().splitlines()
    idx = 0

    nnode, nelem, dim = map(int, lines[idx].split())
    idx += 1

    vertices: List[Tuple[float, float]] = []
    for _ in range(nnode):
        x, y = map(float, lines[idx].split()[:2])
        vertices.append((x, y))
        idx += 1

    nbfgrp = int(lines[idx].split()[0])
    idx += 1

    boundaries: List[BoundaryGroup] = []
    for _ in range(nbfgrp):
        parts = lines[idx].split()
        idx += 1
        nface = int(parts[0])
        nnode_b = int(parts[1])
        title = parts[2] if len(parts) >= 3 else f"BGroup{len(boundaries)+1}"

        faces: List[List[int]] = []
        for _ in range(nface):
            vals = list(map(int, lines[idx].split()))
            if len(vals) != nnode_b:
                raise ValueError(
                    f"Boundary '{title}' expects {nnode_b} nodes/face, got {len(vals)}"
                )
            faces.append(vals)
            idx += 1
        boundaries.append(BoundaryGroup(nface=nface, nnode=nnode_b, title=title, faces=faces))

    element_groups: List[ElementGroup] = []
    nelem_read = 0
    while nelem_read < nelem:
        parts = lines[idx].split()
        idx += 1
        ne_g = int(parts[0])
        p = int(parts[1])
        basis = parts[2]

        if basis != "TriLagrange":
            raise ValueError(f"Unsupported basis '{basis}', only TriLagrange is handled")

        nnode_e = (p + 1) * (p + 2) // 2
        elems: List[List[int]] = []
        for _ in range(ne_g):
            vals = list(map(int, lines[idx].split()))
            if len(vals) != nnode_e:
                raise ValueError(
                    f"Element group p={p} expects {nnode_e} nodes/elem, got {len(vals)}"
                )
            elems.append(vals)
            idx += 1

        nelem_read += ne_g
        element_groups.append(ElementGroup(nelem=ne_g, p=p, basis=basis, elems=elems))

    periodic_groups: List[PeriodicGroup] = []
    if idx < len(lines):
        parts = lines[idx].split()
        if len(parts) >= 2 and parts[1] == "PeriodicGroup":
            npg = int(parts[0])
            idx += 1
            for _ in range(npg):
                h = lines[idx].split()
                idx += 1
                npair = int(h[0])
                kind = h[1] if len(h) > 1 else "Unknown"
                pairs: List[Tuple[int, int]] = []
                for _ in range(npair):
                    a, b = map(int, lines[idx].split()[:2])
                    pairs.append((a, b))
                    idx += 1
                periodic_groups.append(PeriodicGroup(npair=npair, kind=kind, pairs=pairs))

    if idx != len(lines):
        raise ValueError(f"Parsing ended at line {idx+1}, but file has {len(lines)} lines")

    return GriMesh(
        nnode=nnode,
        nelem=nelem,
        dim=dim,
        vertices=vertices,
        boundaries=boundaries,
        element_groups=element_groups,
        periodic_groups=periodic_groups,
    )


def write_gri(mesh: GriMesh, out_path: Path) -> None:
    with out_path.open("w") as f:
        f.write(f"{mesh.nnode} {mesh.nelem} {mesh.dim}\n")
        for x, y in mesh.vertices:
            f.write(f"{x:.15E} {y:.15E}\n")

        f.write(f"{len(mesh.boundaries)}\n")
        for bg in mesh.boundaries:
            f.write(f"{bg.nface} {bg.nnode} {bg.title}\n")
            for face in bg.faces:
                f.write(" ".join(str(v) for v in face) + "\n")

        for eg in mesh.element_groups:
            f.write(f"{eg.nelem} {eg.p} {eg.basis}\n")
            for elem in eg.elems:
                f.write(" ".join(str(v) for v in elem) + "\n")

        if mesh.periodic_groups:
            f.write(f"{len(mesh.periodic_groups)} PeriodicGroup\n")
            for pg in mesh.periodic_groups:
                f.write(f"{pg.npair} {pg.kind}\n")
                for a, b in pg.pairs:
                    f.write(f"{a} {b}\n")


def read_xy_points(path: Path) -> List[Tuple[float, float]]:
    pts: List[Tuple[float, float]] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        x, y = map(float, s.split()[:2])
        pts.append((x, y))
    if len(pts) < 3:
        raise ValueError(f"{path} has too few points ({len(pts)})")
    return pts


class ParametricSpline:
    def __init__(self, pts: List[Tuple[float, float]]):
        arr = np.asarray(pts, dtype=float)
        ds = np.linalg.norm(arr[1:] - arr[:-1], axis=1)
        s = np.concatenate([[0.0], np.cumsum(ds)])
        if s[-1] <= 0:
            raise ValueError("Degenerate curve with zero arclength")
        s /= s[-1]
        self.s = s
        self.csx = CubicSpline(s, arr[:, 0], bc_type="natural")
        self.csy = CubicSpline(s, arr[:, 1], bc_type="natural")

    def xy(self, t: float) -> Tuple[float, float]:
        tt = float(np.clip(t, 0.0, 1.0))
        return float(self.csx(tt)), float(self.csy(tt))

    def project(self, p: Tuple[float, float], nsample: int = 2000) -> Tuple[Tuple[float, float], float]:
        px, py = p
        ts = np.linspace(0.0, 1.0, nsample)
        xs = self.csx(ts)
        ys = self.csy(ts)
        d2 = (xs - px) ** 2 + (ys - py) ** 2
        i0 = int(np.argmin(d2))

        # local bracket around best coarse sample
        i_lo = max(0, i0 - 1)
        i_hi = min(nsample - 1, i0 + 1)
        t_lo = float(ts[i_lo])
        t_hi = float(ts[i_hi])
        if t_hi <= t_lo:
            t_lo, t_hi = 0.0, 1.0

        def obj(t: float) -> float:
            x = float(self.csx(t))
            y = float(self.csy(t))
            return (x - px) * (x - px) + (y - py) * (y - py)

        res = minimize_scalar(obj, bounds=(t_lo, t_hi), method="bounded")
        t_star = float(res.x)
        q = self.xy(t_star)
        return q, float(math.sqrt(obj(t_star)))


def build_edge_group_lookup(mesh: GriMesh) -> Dict[Tuple[int, int], str]:
    lookup: Dict[Tuple[int, int], str] = {}
    for bg in mesh.boundaries:
        if bg.nnode != 2:
            raise ValueError(
                f"Boundary group '{bg.title}' has nnode={bg.nnode}; this script expects edge lists (nnode=2)."
            )
        for n1, n2 in bg.faces:
            a = n1 - 1
            b = n2 - 1
            key = (a, b) if a < b else (b, a)
            lookup[key] = bg.title
    return lookup


def tri_linear_elements(mesh: GriMesh) -> List[List[int]]:
    elems: List[List[int]] = []
    for eg in mesh.element_groups:
        if eg.basis != "TriLagrange" or eg.p != 1:
            raise ValueError(
                f"Input mesh must be linear TriLagrange only. Found group p={eg.p}, basis={eg.basis}."
            )
        for e in eg.elems:
            if len(e) != 3:
                raise ValueError("Linear TriLagrange element must have 3 nodes")
            elems.append(e[:])
    if len(elems) != mesh.nelem:
        raise ValueError(f"Linear element count mismatch: got {len(elems)} vs nelem={mesh.nelem}")
    return elems


def make_curved_mesh_q2(
    mesh: GriMesh,
    lower_curve: ParametricSpline,
    upper_curve: ParametricSpline,
    lower_group: str,
    upper_group: str,
) -> Tuple[GriMesh, Dict[str, float]]:
    edge_to_group = build_edge_group_lookup(mesh)
    linear_elems = tri_linear_elements(mesh)
    verts = list(mesh.vertices)  # 0-based list of coordinates

    # Identify blade boundary edges by group name
    blade_edges: Dict[Tuple[int, int], str] = {}
    for edge, gname in edge_to_group.items():
        if gname == lower_group:
            blade_edges[edge] = "lower"
        elif gname == upper_group:
            blade_edges[edge] = "upper"

    # Curved elements: those touching blade boundary edges
    curved_elem_ids: List[int] = []
    elem_edge_curve_tag: Dict[Tuple[int, int], str] = {}  # (elem_id, local_face_id) -> lower/upper
    for ei, e in enumerate(linear_elems):
        n0, n1, n2 = (e[0] - 1, e[1] - 1, e[2] - 1)  # to 0-based
        face_edges = [
            ((n0, n1) if n0 < n1 else (n1, n0)),
            ((n1, n2) if n1 < n2 else (n2, n1)),
            ((n2, n0) if n2 < n0 else (n0, n2)),
        ]
        touched = False
        for lf, edge in enumerate(face_edges):
            if edge in blade_edges:
                elem_edge_curve_tag[(ei, lf)] = blade_edges[edge]
                touched = True
        if touched:
            curved_elem_ids.append(ei)

    curved_set = set(curved_elem_ids)
    linear_only_ids = [i for i in range(len(linear_elems)) if i not in curved_set]

    # Build global midpoint nodes (shared by all curved elements)
    midpoint_node_id: Dict[Tuple[int, int], int] = {}  # edge -> 1-based node id in output mesh
    proj_dists: List[float] = []

    def get_mid_node(edge: Tuple[int, int], curve_kind: Optional[str]) -> int:
        if edge in midpoint_node_id:
            return midpoint_node_id[edge]

        a, b = edge
        xa, ya = verts[a]
        xb, yb = verts[b]
        mx, my = 0.5 * (xa + xb), 0.5 * (ya + yb)

        if curve_kind == "lower":
            (cx, cy), dist = lower_curve.project((mx, my))
            proj_dists.append(dist)
        elif curve_kind == "upper":
            (cx, cy), dist = upper_curve.project((mx, my))
            proj_dists.append(dist)
        else:
            cx, cy = mx, my
        verts.append((cx, cy))
        out_id = len(verts)  # 1-based node id
        midpoint_node_id[edge] = out_id
        return out_id

    curved_elems_q2: List[List[int]] = []
    for ei in curved_elem_ids:
        e = linear_elems[ei]
        n0, n1, n2 = e  # 1-based vertex ids
        a, b, c = n0 - 1, n1 - 1, n2 - 1

        e01 = (a, b) if a < b else (b, a)  # face 0
        e12 = (b, c) if b < c else (c, b)  # face 1
        e02 = (a, c) if a < c else (c, a)  # face 2

        k01 = elem_edge_curve_tag.get((ei, 0), None)
        k12 = elem_edge_curve_tag.get((ei, 1), None)
        k20 = elem_edge_curve_tag.get((ei, 2), None)

        m01 = get_mid_node(e01, k01)
        m12 = get_mid_node(e12, k12)
        m02 = get_mid_node(e02, k20)

        # p=2 TriLagrange node ordering matching plotgri.m convention:
        # [v0, mid(v0,v1), v1, mid(v0,v2), mid(v1,v2), v2]
        curved_elems_q2.append([n0, m01, n1, m02, m12, n2])

    linear_elems_p1: List[List[int]] = [linear_elems[i] for i in linear_only_ids]

    out_groups: List[ElementGroup] = []
    if curved_elems_q2:
        out_groups.append(
            ElementGroup(
                nelem=len(curved_elems_q2),
                p=2,
                basis="TriLagrange",
                elems=curved_elems_q2,
            )
        )
    if linear_elems_p1:
        out_groups.append(
            ElementGroup(
                nelem=len(linear_elems_p1),
                p=1,
                basis="TriLagrange",
                elems=linear_elems_p1,
            )
        )

    out_mesh = GriMesh(
        nnode=len(verts),
        nelem=mesh.nelem,
        dim=mesh.dim,
        vertices=verts,
        boundaries=mesh.boundaries,
        element_groups=out_groups,
        periodic_groups=mesh.periodic_groups,
    )

    stats = {
        "curved_elements": float(len(curved_elems_q2)),
        "linear_elements": float(len(linear_elems_p1)),
        "new_nodes_added": float(len(verts) - len(mesh.vertices)),
        "curved_boundary_edges_projected": float(len(proj_dists)),
        "projection_dist_mean": float(np.mean(proj_dists) if proj_dists else 0.0),
        "projection_dist_max": float(np.max(proj_dists) if proj_dists else 0.0),
    }
    return out_mesh, stats


def edge_node_map_for_p(p: int) -> List[List[int]]:
    """Return edge node positions (0-based local indices) following plotgri.m logic."""
    if p < 1:
        raise ValueError("p must be >= 1")
    # 1-based build as in plotgri.m
    f1 = list(range(1, p + 2))
    f2: List[int] = []
    v = p + 1
    d = p
    for _ in range(p + 1):
        f2.append(v)
        v += d
        d -= 1
    f3: List[int] = []
    v = 1
    d = p + 1
    for _ in range(p + 1):
        f3.append(v)
        v += d
        d -= 1
    return [[x - 1 for x in f1], [x - 1 for x in f2], [x - 1 for x in f3]]


def sample_quadratic_edge(p0: np.ndarray, pm: np.ndarray, p1: np.ndarray, n: int = 24) -> np.ndarray:
    t = np.linspace(0.0, 1.0, n)
    l0 = 2.0 * (t - 0.5) * (t - 1.0)
    l1 = 4.0 * t * (1.0 - t)
    l2 = 2.0 * t * (t - 0.5)
    xy = (
        l0[:, None] * p0[None, :]
        + l1[:, None] * pm[None, :]
        + l2[:, None] * p1[None, :]
    )
    return xy


def plot_mesh(mesh: GriMesh, title: str, out_png: Path, reference_mesh: Optional[GriMesh] = None) -> None:
    plt = _import_pyplot()
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 5.0), dpi=160)

    # Optional light background mesh (for before/after overlay)
    if reference_mesh is not None:
        for eg in reference_mesh.element_groups:
            if eg.basis != "TriLagrange":
                continue
            edge_map = edge_node_map_for_p(eg.p)
            for elem in eg.elems:
                pts = np.asarray([reference_mesh.vertices[nid - 1] for nid in elem], dtype=float)
                for enodes in edge_map:
                    epts = pts[enodes]
                    ax.plot(epts[:, 0], epts[:, 1], color="#c8c8c8", linewidth=0.45)

    # Current mesh
    for eg in mesh.element_groups:
        if eg.basis != "TriLagrange":
            continue
        edge_map = edge_node_map_for_p(eg.p)
        for elem in eg.elems:
            pts = np.asarray([mesh.vertices[nid - 1] for nid in elem], dtype=float)
            for enodes in edge_map:
                epts = pts[enodes]
                if eg.p == 2 and len(enodes) == 3:
                    curve = sample_quadratic_edge(epts[0], epts[1], epts[2], n=24)
                    ax.plot(curve[:, 0], curve[:, 1], color="k", linewidth=0.55)
                else:
                    ax.plot(epts[:, 0], epts[:, 1], color="k", linewidth=0.55)

    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_before_after(before: GriMesh, after: GriMesh, out_png: Path) -> None:
    plt = _import_pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=160, sharex=True, sharey=True)

    for ax, mesh, title in [
        (axes[0], before, "Before (linear)"),
        (axes[1], after, "After (q=2 curved near blade)"),
    ]:
        for eg in mesh.element_groups:
            if eg.basis != "TriLagrange":
                continue
            edge_map = edge_node_map_for_p(eg.p)
            for elem in eg.elems:
                pts = np.asarray([mesh.vertices[nid - 1] for nid in elem], dtype=float)
                for enodes in edge_map:
                    epts = pts[enodes]
                    if eg.p == 2 and len(enodes) == 3:
                        curve = sample_quadratic_edge(epts[0], epts[1], epts[2], n=24)
                        ax.plot(curve[:, 0], curve[:, 1], color="k", linewidth=0.45)
                    else:
                        ax.plot(epts[:, 0], epts[:, 1], color="k", linewidth=0.45)
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.grid(False)

    axes[0].set_ylabel("y")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> None:
    # Manual configuration: edit these paths/names directly.
    base_dir = Path(__file__).resolve().parent
    input_mesh = base_dir / "mesh3.gri"
    lower_blade_file = base_dir / "bladelower.txt"
    upper_blade_file = base_dir / "bladeupper.txt"
    output_mesh = base_dir / "mesh3_curved.gri"
    lower_group = "BGroup2"
    upper_group = "BGroup6"
    upper_shift_y = -18.0

    if not input_mesh.exists():
        raise FileNotFoundError(f"Input mesh not found: {input_mesh}")
    if not lower_blade_file.exists():
        raise FileNotFoundError(f"Lower blade file not found: {lower_blade_file}")
    if not upper_blade_file.exists():
        raise FileNotFoundError(f"Upper blade file not found: {upper_blade_file}")

    in_mesh = parse_gri(input_mesh)
    lower_pts = read_xy_points(lower_blade_file)
    upper_pts = read_xy_points(upper_blade_file)

    # Bring upper blade into this mesh's periodic copy by default.
    upper_pts_shifted = [(x, y + upper_shift_y) for (x, y) in upper_pts]

    lower_curve = ParametricSpline(lower_pts)
    upper_curve = ParametricSpline(upper_pts_shifted)

    out_mesh, stats = make_curved_mesh_q2(
        in_mesh,
        lower_curve=lower_curve,
        upper_curve=upper_curve,
        lower_group=lower_group,
        upper_group=upper_group,
    )

    write_gri(out_mesh, output_mesh)

    print("Curved mesh generation complete.")
    print(f"  input mesh      : {input_mesh}")
    print(f"  output mesh     : {output_mesh}")
    print("  stats:")
    for k, v in stats.items():
        if abs(v - int(v)) < 1e-12:
            print(f"    - {k}: {int(v)}")
        else:
            print(f"    - {k}: {v:.6e}")


if __name__ == "__main__":
    main()
