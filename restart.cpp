#include "restart.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <utility>

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
