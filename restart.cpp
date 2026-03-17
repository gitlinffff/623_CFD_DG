#include "restart.hpp"
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
            mkdir(sub.c_str(), 0755);
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
            nodes.push_back({0.0, 0.0});
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

struct LinearTriMap {
    double x1;
    double y1;
    double xr;
    double yr;
    double xs;
    double ys;
    double det;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

bool build_linear_tri_map(const GriMesh& mesh, int elem, LinearTriMap& map) {
    const int* tri = &mesh.E[elem * 3];
    const int n1 = tri[0];
    const int n2 = tri[1];
    const int n3 = tri[2];

    const double x1 = mesh.V[n1 * 2 + 0];
    const double y1 = mesh.V[n1 * 2 + 1];
    const double x2 = mesh.V[n2 * 2 + 0];
    const double y2 = mesh.V[n2 * 2 + 1];
    const double x3 = mesh.V[n3 * 2 + 0];
    const double y3 = mesh.V[n3 * 2 + 1];

    map.x1 = x1;
    map.y1 = y1;
    map.xr = x2 - x1;
    map.yr = y2 - y1;
    map.xs = x3 - x1;
    map.ys = y3 - y1;
    map.det = map.xr * map.ys - map.yr * map.xs;
    if (std::fabs(map.det) < 1e-14) return false;

    map.xmin = std::min(x1, std::min(x2, x3));
    map.xmax = std::max(x1, std::max(x2, x3));
    map.ymin = std::min(y1, std::min(y2, y3));
    map.ymax = std::max(y1, std::max(y2, y3));
    return true;
}

inline void ref_to_phys(const LinearTriMap& map, double r, double s, double& x, double& y) {
    x = map.x1 + map.xr * r + map.xs * s;
    y = map.y1 + map.yr * r + map.ys * s;
}

inline void phys_to_ref(const LinearTriMap& map, double x, double y, double& r, double& s) {
    const double dx = x - map.x1;
    const double dy = y - map.y1;
    r = (dx * map.ys - dy * map.xs) / map.det;
    s = (map.xr * dy - map.yr * dx) / map.det;
}

inline bool inside_ref_triangle(double r, double s, double tol) {
    return (r >= -tol) && (s >= -tol) && (r + s <= 1.0 + tol);
}

int find_source_elem_for_point(const std::vector<LinearTriMap>& src_maps, double x, double y) {
    const double tol = 1e-10;
    for (int e = 0; e < (int)src_maps.size(); ++e) {
        const LinearTriMap& m = src_maps[e];
        if (x < m.xmin - tol || x > m.xmax + tol || y < m.ymin - tol || y > m.ymax + tol)
            continue;
        double r = 0.0, s = 0.0;
        phys_to_ref(m, x, y, r, s);
        if (inside_ref_triangle(r, s, 1e-8)) return e;
    }
    return -1;
}

}

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
            std::getline(in, tail);
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
    if (target_order != src.order) {
        err_msg = "cross-mesh restart requires the same DG order.";
        return false;
    }
    if (src_mesh.Ne != src.Ne) {
        err_msg = "source mesh element count does not match restart Ne.";
        return false;
    }

    const int np = num_basis(target_order);
    if (src.Np != np) {
        err_msg = "restart Np is inconsistent with restart order.";
        return false;
    }
    if ((int)src.U.size() != 4 * src.Ne * src.Np) {
        err_msg = "restart data size is inconsistent.";
        return false;
    }

    std::vector<std::pair<double, double>> ref_nodes;
    if (!reference_nodes(target_order, ref_nodes) || (int)ref_nodes.size() != np) {
        err_msg = "failed to build reference nodes for cross-mesh restart.";
        return false;
    }

    std::vector<LinearTriMap> src_maps(src_mesh.Ne);
    for (int e = 0; e < src_mesh.Ne; ++e) {
        if (!build_linear_tri_map(src_mesh, e, src_maps[e])) {
            err_msg = "degenerate source triangle found while building cross-mesh restart.";
            return false;
        }
    }

    std::vector<LinearTriMap> dst_maps(target_mesh.Ne);
    for (int e = 0; e < target_mesh.Ne; ++e) {
        if (!build_linear_tri_map(target_mesh, e, dst_maps[e])) {
            err_msg = "degenerate target triangle found while building cross-mesh restart.";
            return false;
        }
    }

    std::vector<int> parent(target_mesh.Ne, -1);
    for (int e = 0; e < target_mesh.Ne; ++e) {
        double xc = 0.0, yc = 0.0;
        ref_to_phys(dst_maps[e], 1.0 / 3.0, 1.0 / 3.0, xc, yc);
        parent[e] = find_source_elem_for_point(src_maps, xc, yc);
        if (parent[e] < 0) {
            std::ostringstream oss;
            oss << "failed to find source element for target element " << e
                << " (centroid at x=" << xc << ", y=" << yc << ").";
            err_msg = oss.str();
            return false;
        }
    }

    double* phi = nullptr;
    for (int e = 0; e < target_mesh.Ne; ++e) {
        const int src_elem = parent[e];
        const LinearTriMap& src_map = src_maps[src_elem];
        const LinearTriMap& dst_map = dst_maps[e];
        for (int m = 0; m < np; ++m) {
            const double rt = ref_nodes[m].first;
            const double st = ref_nodes[m].second;
            double x = 0.0, y = 0.0, rs = 0.0, ss = 0.0;
            ref_to_phys(dst_map, rt, st, x, y);
            phys_to_ref(src_map, x, y, rs, ss);

            double xref[2] = {rs, ss};
            shapeL(xref, src.order, &phi);
            for (int var = 0; var < 4; ++var) {
                const int src_base = (var * src_mesh.Ne + src_elem) * np;
                const int dst_idx = (var * target_mesh.Ne + e) * np + m;
                double val = 0.0;
                for (int j = 0; j < np; ++j) {
                    val += phi[j] * src.U[src_base + j];
                }
                U_target[dst_idx] = val;
            }
        }
    }
    if (phi) free(phi);

    return true;
}
