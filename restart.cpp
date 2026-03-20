#include "restart.hpp"
#include "geometry.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <utility>
#include <vector>

extern "C" {
    void shapeL(double* xref, int p, double** pphi);
}

namespace {

int num_basis(int order) {
    return (order + 1) * (order + 2) / 2;
}

void make_dirs(const std::string& path) {
    for (size_t i = 1; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            std::string sub = path.substr(0, i);
            mkdir(sub.c_str(), 0755); // ignore EEXIST
        }
    }
}

std::string ltrim(const std::string& s) {
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r')) ++i;
    return s.substr(i);
}

bool reference_nodes(int order, std::vector<std::pair<double, double>>& nodes) {
    nodes.clear();
    switch (order) {
        case 0:
            // For p=0 we sample at centroid to avoid boundary-point ambiguity
            // during cross-mesh point location.
            nodes.push_back({1.0 / 3.0, 1.0 / 3.0});
            return true;
        case 1:
            nodes.push_back({0.0, 0.0});
            nodes.push_back({1.0, 0.0});
            nodes.push_back({0.0, 1.0});
            return true;
        case 2:
            nodes.push_back({0.0, 0.0});
            nodes.push_back({1.0, 0.0});
            nodes.push_back({0.0, 1.0});
            nodes.push_back({0.5, 0.5});
            nodes.push_back({0.0, 0.5});
            nodes.push_back({0.5, 0.0});
            return true;
        case 3:
            nodes.push_back({0.0, 0.0});
            nodes.push_back({1.0, 0.0});
            nodes.push_back({0.0, 1.0});
            nodes.push_back({2.0 / 3.0, 1.0 / 3.0});
            nodes.push_back({1.0 / 3.0, 2.0 / 3.0});
            nodes.push_back({0.0, 2.0 / 3.0});
            nodes.push_back({0.0, 1.0 / 3.0});
            nodes.push_back({1.0 / 3.0, 0.0});
            nodes.push_back({2.0 / 3.0, 0.0});
            nodes.push_back({1.0 / 3.0, 1.0 / 3.0});
            return true;
        default:
            return false;
    }
}

bool interpolation_matrix(int src_order, int dst_order, std::vector<double>& T) {
    std::vector<std::pair<double, double>> dst_nodes;
    if (!reference_nodes(dst_order, dst_nodes)) return false;

    const int np_src = num_basis(src_order);
    const int np_dst = num_basis(dst_order);
    if ((int)dst_nodes.size() != np_dst) return false;

    T.assign(np_dst * np_src, 0.0);
    double* phi = nullptr;
    for (int m = 0; m < np_dst; ++m) {
        double xref[2] = {dst_nodes[m].first, dst_nodes[m].second};
        shapeL(xref, src_order, &phi);
        for (int j = 0; j < np_src; ++j) {
            T[m * np_src + j] = phi[j];
        }
    }
    if (phi) free(phi);
    return true;
}

inline bool inside_ref_triangle(double r, double s, double tol) {
    return (r >= -tol) && (s >= -tol) && (r + s <= 1.0 + tol);
}

struct BBox {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

void get_elem_geom_nodes(const GriMesh& mesh, int elem, const int*& conn, int& nnode) {
    if ((int)mesh.ElemGeomOrder.size() == mesh.Ne &&
        (int)mesh.ElemGeomOffset.size() == mesh.Ne + 1) {
        int q = mesh.ElemGeomOrder[elem];
        int nn = geom_num_nodes(q);
        int off0 = mesh.ElemGeomOffset[elem];
        int off1 = mesh.ElemGeomOffset[elem + 1];
        if (nn > 0 && off0 >= 0 && off1 <= (int)mesh.ElemGeomConn.size() && off1 - off0 == nn) {
            conn = &mesh.ElemGeomConn[off0];
            nnode = nn;
            return;
        }
    }
    conn = &mesh.E[elem * 3];
    nnode = 3;
}

void build_elem_bbox(const GriMesh& mesh, int elem, BBox& box) {
    const int* conn = nullptr;
    int nnode = 0;
    get_elem_geom_nodes(mesh, elem, conn, nnode);

    double xmin = 1e300, xmax = -1e300;
    double ymin = 1e300, ymax = -1e300;
    for (int i = 0; i < nnode; ++i) {
        int nid = conn[i];
        double x = mesh.V[nid * 2 + 0];
        double y = mesh.V[nid * 2 + 1];
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
    }

    const double eps = 1e-12;
    box.xmin = xmin - eps;
    box.xmax = xmax + eps;
    box.ymin = ymin - eps;
    box.ymax = ymax + eps;
}

bool point_in_bbox(const BBox& b, double x, double y) {
    return (x >= b.xmin && x <= b.xmax && y >= b.ymin && y <= b.ymax);
}

int find_source_elem_for_point(const GriMesh& src_mesh,
                               const std::vector<BBox>& src_bbox,
                               double x,
                               double y,
                               int preferred_elem = -1) {
    const double tol = 1e-8;

    if (preferred_elem >= 0 && preferred_elem < src_mesh.Ne) {
        if (point_in_bbox(src_bbox[(size_t)preferred_elem], x, y)) {
            double r = 0.0, s = 0.0;
            if (phys_to_ref(src_mesh, preferred_elem, x, y, r, s) && inside_ref_triangle(r, s, tol))
                return preferred_elem;
        }
    }

    for (int e = 0; e < src_mesh.Ne; ++e) {
        if (!point_in_bbox(src_bbox[(size_t)e], x, y)) continue;
        double r = 0.0, s = 0.0;
        if (!phys_to_ref(src_mesh, e, x, y, r, s)) continue;
        if (inside_ref_triangle(r, s, tol)) return e;
    }
    return -1;
}

} // namespace

bool write_restart_dat(const std::string& filepath,
                       const GriMesh& mesh,
                       int order,
                       const double* U,
                       const std::string& mesh_file) {
    const int Np = num_basis(order);
    const long long nTot = 4LL * mesh.Ne * Np;
    if (nTot <= 0) return false;

    const size_t slash = filepath.find_last_of('/');
    if (slash != std::string::npos) {
        make_dirs(filepath.substr(0, slash));
    }

    std::ofstream out(filepath.c_str());
    if (!out) return false;

    out << "DG_RESTART_V1\n";
    out << "order " << order << "\n";
    out << "Ne " << mesh.Ne << "\n";
    out << "Np " << Np << "\n";
    out << "nvars 4\n";
    out << "mesh_file " << mesh_file << "\n";
    out << "nTot " << nTot << "\n";
    out << std::setprecision(17);
    for (long long i = 0; i < nTot; ++i) {
        out << U[i] << "\n";
    }
    return true;
}

bool read_restart_dat(const std::string& filepath, RestartData& data) {
    std::ifstream in(filepath.c_str());
    if (!in) return false;

    data = RestartData();

    std::string magic;
    in >> magic;
    if (!in || magic != "DG_RESTART_V1") return false;

    std::string key;
    long long nTot = -1;
    while (in >> key) {
        if (key == "order") {
            in >> data.order;
        } else if (key == "Ne") {
            in >> data.Ne;
        } else if (key == "Np") {
            in >> data.Np;
        } else if (key == "nvars") {
            in >> data.nvars;
        } else if (key == "mesh_file") {
            std::string tail;
            std::getline(in, tail);
            data.mesh_file = ltrim(tail);
        } else if (key == "nTot") {
            in >> nTot;
            break;
        } else {
            std::string tail;
            std::getline(in, tail); // skip unknown line
        }
    }

    if (!in || data.order < 0 || data.Ne <= 0 || data.Np <= 0 || data.nvars != 4 || nTot <= 0)
        return false;

    const int expectedNp = num_basis(data.order);
    if (expectedNp != data.Np) return false;
    if (nTot != 4LL * data.Ne * data.Np) return false;

    data.U.resize((size_t)nTot);
    for (long long i = 0; i < nTot; ++i) {
        in >> data.U[(size_t)i];
        if (!in) return false;
    }
    return true;
}

bool initialize_from_restart(const RestartData& src,
                             int target_order,
                             int target_Ne,
                             double* U_target,
                             std::string& err_msg) {
    err_msg.clear();
    if (src.nvars != 4) {
        err_msg = "restart nvars must be 4.";
        return false;
    }
    if (src.Ne != target_Ne) {
        err_msg = "restart Ne does not match current mesh Ne.";
        return false;
    }
    const int np_src = num_basis(src.order);
    const int np_dst = num_basis(target_order);
    if (np_src != src.Np) {
        err_msg = "restart Np is inconsistent with restart order.";
        return false;
    }
    if ((int)src.U.size() != 4 * src.Ne * src.Np) {
        err_msg = "restart data size is inconsistent.";
        return false;
    }

    std::vector<double> T;
    if (!interpolation_matrix(src.order, target_order, T)) {
        err_msg = "failed to build interpolation matrix for restart prolongation.";
        return false;
    }

    for (int var = 0; var < 4; ++var) {
        for (int k = 0; k < target_Ne; ++k) {
            const int src_base = (var * target_Ne + k) * np_src;
            const int dst_base = (var * target_Ne + k) * np_dst;
            for (int m = 0; m < np_dst; ++m) {
                double val = 0.0;
                const int row = m * np_src;
                for (int j = 0; j < np_src; ++j) {
                    val += T[row + j] * src.U[src_base + j];
                }
                U_target[dst_base + m] = val;
            }
        }
    }

    return true;
}

bool initialize_from_restart_cross_mesh(const RestartData& src,
                                        const GriMesh& src_mesh,
                                        int target_order,
                                        const GriMesh& target_mesh,
                                        double* U_target,
                                        std::string& err_msg) {
    err_msg.clear();

    if (src.nvars != 4) {
        err_msg = "restart nvars must be 4.";
        return false;
    }
    if (src_mesh.Ne != src.Ne) {
        err_msg = "source mesh element count does not match restart Ne.";
        return false;
    }

    const int np_src = num_basis(src.order);
    const int np_dst = num_basis(target_order);
    if (src.Np != np_src) {
        err_msg = "restart Np is inconsistent with restart order.";
        return false;
    }
    if ((int)src.U.size() != 4 * src.Ne * src.Np) {
        err_msg = "restart data size is inconsistent.";
        return false;
    }

    std::vector<std::pair<double, double>> dst_ref_nodes;
    if (!reference_nodes(target_order, dst_ref_nodes) || (int)dst_ref_nodes.size() != np_dst) {
        err_msg = "failed to build target reference nodes for cross-mesh restart.";
        return false;
    }

    std::vector<BBox> src_bbox((size_t)src_mesh.Ne);
    for (int e = 0; e < src_mesh.Ne; ++e) {
        build_elem_bbox(src_mesh, e, src_bbox[(size_t)e]);
    }

    // Determine source parent by target centroid.
    std::vector<int> parent((size_t)target_mesh.Ne, -1);
    for (int e = 0; e < target_mesh.Ne; ++e) {
        GeomEval g;
        eval_geometry(target_mesh, e, 1.0 / 3.0, 1.0 / 3.0, g);
        int src_elem = find_source_elem_for_point(src_mesh, src_bbox, g.x, g.y, -1);
        if (src_elem < 0) {
            std::ostringstream oss;
            oss << "failed to locate parent source element for target element " << e
                << " (centroid x=" << g.x << ", y=" << g.y << ").";
            err_msg = oss.str();
            return false;
        }
        parent[(size_t)e] = src_elem;
    }

    double* phi = nullptr;
    for (int e = 0; e < target_mesh.Ne; ++e) {
        int src_elem_hint = parent[(size_t)e];
        for (int m = 0; m < np_dst; ++m) {
            double rt = dst_ref_nodes[(size_t)m].first;
            double st = dst_ref_nodes[(size_t)m].second;

            GeomEval gt;
            eval_geometry(target_mesh, e, rt, st, gt);

            // Try parent hint first, fallback to global search if needed.
            int src_elem = src_elem_hint;
            double rs = 0.0, ss = 0.0;
            bool ok = phys_to_ref(src_mesh, src_elem, gt.x, gt.y, rs, ss) &&
                      inside_ref_triangle(rs, ss, 1e-8);
            if (!ok) {
                src_elem = find_source_elem_for_point(src_mesh, src_bbox, gt.x, gt.y, src_elem_hint);
                if (src_elem < 0) {
                    std::ostringstream oss;
                    oss << "failed to locate source element for target interpolation point"
                        << " (target elem " << e << ", node " << m
                        << ", x=" << gt.x << ", y=" << gt.y << ").";
                    err_msg = oss.str();
                    if (phi) free(phi);
                    return false;
                }
                if (!phys_to_ref(src_mesh, src_elem, gt.x, gt.y, rs, ss)) {
                    err_msg = "phys_to_ref failed during cross-mesh restart interpolation.";
                    if (phi) free(phi);
                    return false;
                }
            }

            double xref_src[2] = {rs, ss};
            shapeL(xref_src, src.order, &phi);

            for (int var = 0; var < 4; ++var) {
                const int src_base = (var * src_mesh.Ne + src_elem) * np_src;
                const int dst_idx = (var * target_mesh.Ne + e) * np_dst + m;
                double val = 0.0;
                for (int j = 0; j < np_src; ++j) {
                    val += phi[j] * src.U[src_base + j];
                }
                U_target[dst_idx] = val;
            }
        }
    }

    if (phi) free(phi);
    return true;
}
