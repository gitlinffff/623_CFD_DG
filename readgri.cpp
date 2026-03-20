#include "readgri.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <cmath>
#include <cctype>
#include <cstring>
#include <cstdlib>

namespace {

using EdgeKey = std::pair<int, int>;

struct VisitedEntry {
    int elemL;   /* 0-based */
    int faceL;   /* 0-based local face */
};

bool read_line_ints(std::ifstream& f, int* out, int n) {
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    for (int i = 0; i < n; ++i)
        if (!(ss >> out[i])) return false;
    return true;
}

bool read_line_doubles(std::ifstream& f, double* out, int n) {
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    for (int i = 0; i < n; ++i)
        if (!(ss >> out[i])) return false;
    return true;
}

bool read_line_var_ints(std::ifstream& f, int n, std::vector<int>& out) {
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    out.resize(n);
    for (int i = 0; i < n; ++i) {
        if (!(ss >> out[i])) return false;
    }
    return true;
}

double norm2(double x, double y) {
    return std::sqrt(x * x + y * y);
}

double signed_twice_area(const GriMesh& mesh, int n0, int n1, int n2) {
    const double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
    const double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
    const double x2 = mesh.V[n2 * 2 + 0], y2 = mesh.V[n2 * 2 + 1];
    return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
}

void face_outward_normal(const GriMesh& mesh, int elem, int face,
                         double& nx, double& ny, double& len) {
    const int* tri = &mesh.E[elem * 3];
    int n0 = tri[face];
    int n1 = tri[(face + 1) % 3];
    int n2 = tri[(face + 2) % 3];

    double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
    double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
    double x2 = mesh.V[n2 * 2 + 0], y2 = mesh.V[n2 * 2 + 1];

    double ex = x1 - x0;
    double ey = y1 - y0;
    len = norm2(ex, ey);
    if (len < 1e-20) len = 1e-20;
    ex /= len;
    ey /= len;

    nx = ey;
    ny = -ex;

    // Enforce outward orientation with respect to element centroid.
    double mx = 0.5 * (x0 + x1);
    double my = 0.5 * (y0 + y1);
    double cx = (x0 + x1 + x2) / 3.0;
    double cy = (y0 + y1 + y2) / 3.0;
    double dot = nx * (mx - cx) + ny * (my - cy);
    if (dot < 0.0) {
        nx = -nx;
        ny = -ny;
    }
}

void reorder_p2_plotgri_to_shape(const std::vector<int>& in_plot, std::vector<int>& out_shape) {
    // Input (.gri / plotgri order): [v0, m01, v1, m02, m12, v2]
    // shapeL p=2 order:            [v0, v1,  v2, m12, m20, m01]
    out_shape.resize(6);
    out_shape[0] = in_plot[0]; // v0
    out_shape[1] = in_plot[2]; // v1
    out_shape[2] = in_plot[5]; // v2
    out_shape[3] = in_plot[4]; // m12
    out_shape[4] = in_plot[3]; // m20
    out_shape[5] = in_plot[1]; // m01
}

void enforce_ccw_local_ordering(const GriMesh& mesh, int p, std::vector<int>& shape_conn) {
    if (shape_conn.size() < 3) return;
    if (signed_twice_area(mesh, shape_conn[0], shape_conn[1], shape_conn[2]) >= 0.0) return;

    if (p == 1) {
        // [v0, v1, v2] -> [v0, v2, v1]
        std::swap(shape_conn[1], shape_conn[2]);
        return;
    }

    if (p == 2 && shape_conn.size() == 6) {
        // shapeL p=2 order: [v0, v1, v2, m12, m20, m01]
        // Flipped order:     [v0, v2, v1, m12, m01, m20]
        std::vector<int> t(6);
        t[0] = shape_conn[0];
        t[1] = shape_conn[2];
        t[2] = shape_conn[1];
        t[3] = shape_conn[3];
        t[4] = shape_conn[5];
        t[5] = shape_conn[4];
        shape_conn.swap(t);
        return;
    }
}

}  // namespace

bool read_gri(const char* fname, GriMesh& mesh) {
    std::ifstream f(fname);
    if (!f.is_open()) return false;

    {
        int buf[3];
        if (!read_line_ints(f, buf, 3)) return false;
        mesh.Nn = buf[0];
        mesh.Ne = buf[1];
    }

    mesh.V.resize(mesh.Nn * 2);
    for (int n = 0; n < mesh.Nn; ++n) {
        if (!read_line_doubles(f, &mesh.V[n * 2], 2)) return false;
    }

    int NB;
    {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        if (!(ss >> NB)) return false;
    }

    std::map<EdgeKey, int> BgroupOfEdge;
    mesh.Bname.clear();
    mesh.Bname.reserve((size_t)NB);
    for (int i = 0; i < NB; ++i) {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        int Nb = 0, nnode_b = 0;
        std::string name;
        if (!(ss >> Nb >> nnode_b >> name)) return false;
        if (nnode_b < 2) return false;

        mesh.Bname.push_back(name);
        for (int n = 0; n < Nb; ++n) {
            std::vector<int> bnodes;
            if (!read_line_var_ints(f, nnode_b, bnodes)) return false;
            for (int& id : bnodes) --id;
            int n0 = bnodes.front();
            int n1 = bnodes.back();
            int a = n0 < n1 ? n0 : n1;
            int b = n0 < n1 ? n1 : n0;
            BgroupOfEdge[EdgeKey(a, b)] = i + 1; // bgroup is 1-based
        }
    }

    mesh.E.resize(mesh.Ne * 3);
    mesh.ElemGeomOrder.assign(mesh.Ne, 1);
    mesh.ElemGeomOffset.assign(mesh.Ne + 1, 0);
    std::vector<int> geom_conn;
    geom_conn.reserve((size_t)mesh.Ne * 3);

    int elem_cursor = 0;
    const char* disable_ccw_fix = std::getenv("DG_DISABLE_CCW_FIX");
    const bool do_ccw_fix = !(disable_ccw_fix && std::strcmp(disable_ccw_fix, "1") == 0);

    while (elem_cursor < mesh.Ne) {
        std::string line;
        if (!std::getline(f, line)) return false;
        if (line.empty()) continue;

        std::istringstream ss(line);
        int ne_g = 0, p = 0;
        std::string basis;
        if (!(ss >> ne_g >> p >> basis)) return false;
        if (ne_g <= 0) return false;
        if (basis != "TriLagrange") return false;
        if (p < 1 || p > 2) return false; // current solver supports geometric q=1/2

        const int nnode_e = (p + 1) * (p + 2) / 2;
        for (int n = 0; n < ne_g; ++n) {
            if (elem_cursor >= mesh.Ne) return false;

            std::vector<int> raw_conn;
            if (!read_line_var_ints(f, nnode_e, raw_conn)) return false;
            for (int& id : raw_conn) --id;

            std::vector<int> shape_conn;
            if (p == 1) {
                if ((int)raw_conn.size() != 3) return false;
                shape_conn = raw_conn;
            } else { // p == 2
                if ((int)raw_conn.size() != 6) return false;
                reorder_p2_plotgri_to_shape(raw_conn, shape_conn);
            }

            if (do_ccw_fix)
                enforce_ccw_local_ordering(mesh, p, shape_conn);

            mesh.E[elem_cursor * 3 + 0] = shape_conn[0];
            mesh.E[elem_cursor * 3 + 1] = shape_conn[1];
            mesh.E[elem_cursor * 3 + 2] = shape_conn[2];

            mesh.ElemGeomOrder[elem_cursor] = p;
            mesh.ElemGeomOffset[elem_cursor] = (int)geom_conn.size();
            geom_conn.insert(geom_conn.end(), shape_conn.begin(), shape_conn.end());
            ++elem_cursor;
        }
    }
    if (elem_cursor != mesh.Ne) return false;
    mesh.ElemGeomOffset[mesh.Ne] = (int)geom_conn.size();
    mesh.ElemGeomConn = std::move(geom_conn);

    std::vector<std::map<int, int>> PG;
    {
        std::string line;
        if (std::getline(f, line)) {
            bool has_non_ws = false;
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) {
                    has_non_ws = true;
                    break;
                }
            }
            if (has_non_ws) {
                std::istringstream ss(line);
                int NPG = 0;
                if (!(ss >> NPG)) return false;
                if (NPG < 0) return false;
                PG.resize((size_t)NPG);
                for (int i = 0; i < NPG; ++i) {
                    std::string hline;
                    if (!std::getline(f, hline)) return false;
                    std::istringstream hs(hline);
                    int Npgn = 0;
                    std::string type;
                    if (!(hs >> Npgn >> type)) return false;
                    if (Npgn < 0) return false;
                    for (int n = 0; n < Npgn; ++n) {
                        int pair_buf[2];
                        if (!read_line_ints(f, pair_buf, 2)) return false;
                        int n1 = pair_buf[0] - 1;
                        int n2 = pair_buf[1] - 1;
                        PG[(size_t)i][n1] = n2;
                        PG[(size_t)i][n2] = n1;
                    }
                }
            }
        }
    }
    f.close();

    mesh.Area.resize(mesh.Ne);
    std::vector<int> I2E;
    std::vector<double> In, In_len;
    std::vector<int> B2E;
    std::vector<double> Bn, Bn_len;
    std::map<EdgeKey, VisitedEntry> visited;

    for (int i = 0; i < mesh.Ne; ++i) {
        int* tri = &mesh.E[i * 3];
        double x1 = mesh.V[tri[0] * 2 + 0], y1 = mesh.V[tri[0] * 2 + 1];
        double x2 = mesh.V[tri[1] * 2 + 0], y2 = mesh.V[tri[1] * 2 + 1];
        double x3 = mesh.V[tri[2] * 2 + 0], y3 = mesh.V[tri[2] * 2 + 1];
        mesh.Area[i] = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));

        for (int j = 0; j < 3; ++j) {
            int na = tri[j], nb = tri[(j + 1) % 3];
            int end0 = na < nb ? na : nb;
            int end1 = na < nb ? nb : na;
            EdgeKey faceID(end0, end1);

            auto itB = BgroupOfEdge.find(faceID);
            if (itB != BgroupOfEdge.end()) {
                bool periodic = false;
                int pg_other0 = -1, pg_other1 = -1;
                for (const auto& m : PG) {
                    auto it0 = m.find(end0);
                    auto it1 = m.find(end1);
                    if (it0 != m.end() && it1 != m.end()) {
                        periodic = true;
                        pg_other0 = it0->second;
                        pg_other1 = it1->second;
                        break;
                    }
                }

                if (periodic) {
                    int o0 = pg_other0 < pg_other1 ? pg_other0 : pg_other1;
                    int o1 = pg_other0 < pg_other1 ? pg_other1 : pg_other0;
                    EdgeKey otherID(o0, o1);
                    auto itV = visited.find(otherID);
                    if (itV != visited.end()) {
                        VisitedEntry& ve = itV->second;
                        I2E.push_back(ve.elemL);
                        I2E.push_back(ve.faceL);
                        I2E.push_back(i);
                        I2E.push_back(j);

                        double nx, ny, len;
                        face_outward_normal(mesh, ve.elemL, ve.faceL, nx, ny, len);
                        In.push_back(nx);
                        In.push_back(ny);
                        In_len.push_back(len);
                        visited.erase(itV);
                    } else {
                        visited[faceID] = { i, j };
                    }
                } else {
                    int bgroup = itB->second;
                    B2E.push_back(i);
                    B2E.push_back(j);
                    B2E.push_back(bgroup);

                    double nx, ny, len;
                    face_outward_normal(mesh, i, j, nx, ny, len);
                    Bn.push_back(nx);
                    Bn.push_back(ny);
                    Bn_len.push_back(len);
                }
            } else {
                auto itV = visited.find(faceID);
                if (itV != visited.end()) {
                    VisitedEntry& ve = itV->second;
                    I2E.push_back(ve.elemL);
                    I2E.push_back(ve.faceL);
                    I2E.push_back(i);
                    I2E.push_back(j);

                    double nx, ny, len;
                    face_outward_normal(mesh, ve.elemL, ve.faceL, nx, ny, len);
                    In.push_back(nx);
                    In.push_back(ny);
                    In_len.push_back(len);
                    visited.erase(itV);
                } else {
                    visited[faceID] = { i, j };
                }
            }
        }
    }

    mesh.num_interior_faces = (int)I2E.size() / 4;
    mesh.I2E = std::move(I2E);
    mesh.In = std::move(In);
    mesh.In_len = std::move(In_len);
    mesh.num_boundary_faces = (int)B2E.size() / 3;
    mesh.B2E = std::move(B2E);
    mesh.Bn = std::move(Bn);
    mesh.Bn_len = std::move(Bn_len);

    return true;
}
