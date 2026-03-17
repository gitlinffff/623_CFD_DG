#include "readgri.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <cmath>
#include <cstring>
#include <cstdlib>

namespace {

using EdgeKey = std::pair<int, int>;

struct VisitedEntry {
    int elemL;
    int faceL;
};

struct BoundaryEdgeInfo {
    int bgroup;
    std::vector<int> nodes;
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

double norm2(double x, double y) {
    return std::sqrt(x * x + y * y);
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

}

bool read_gri(const char* fname, GriMesh& mesh) {
    std::ifstream f(fname);
    if (!f.is_open()) return false;

    int dim;
    {
        int buf[3];
        if (!read_line_ints(f, buf, 3)) return false;
        mesh.Nn = buf[0];
        mesh.Ne = buf[1];
        dim = buf[2];
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

    std::map<EdgeKey, BoundaryEdgeInfo> B;
    mesh.Bname.clear();
    for (int i = 0; i < NB; ++i) {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        int Nb, nnode;
        std::string name;
        if (!(ss >> Nb >> nnode >> name)) return false;
        if (nnode < 2) return false;
        mesh.Bname.push_back(name);
        int bgroup = (int)mesh.Bname.size();
        for (int n = 0; n < Nb; ++n) {
            std::string line2;
            if (!std::getline(f, line2)) return false;
            std::istringstream ss2(line2);
            std::vector<int> nodes((size_t)nnode);
            for (int k = 0; k < nnode; ++k) {
                if (!(ss2 >> nodes[(size_t)k])) return false;
                nodes[(size_t)k]--;
            }
            int n1 = nodes[0];
            int n2 = nodes[1];
            int a = n1 < n2 ? n1 : n2;
            int b = n1 < n2 ? n2 : n1;
            BoundaryEdgeInfo info;
            info.bgroup = bgroup;
            info.nodes = std::move(nodes);
            B[EdgeKey(a, b)] = std::move(info);
        }
    }

    mesh.E.resize(mesh.Ne * 3);
    int Ne0 = 0;
    while (Ne0 < mesh.Ne) {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        int ne;
        if (!(ss >> ne)) return false;
        for (int n = 0; n < ne; ++n) {
            if (!read_line_ints(f, &mesh.E[(Ne0 + n) * 3], 3)) return false;
            for (int k = 0; k < 3; ++k)
                mesh.E[(Ne0 + n) * 3 + k]--;
        }
        Ne0 += ne;
    }

    const char* disable_ccw_fix = std::getenv("DG_DISABLE_CCW_FIX");
    if (!(disable_ccw_fix && std::strcmp(disable_ccw_fix, "1") == 0)) {
        for (int e = 0; e < mesh.Ne; ++e) {
            int* tri = &mesh.E[e * 3];
            double x1 = mesh.V[tri[0] * 2 + 0], y1 = mesh.V[tri[0] * 2 + 1];
            double x2 = mesh.V[tri[1] * 2 + 0], y2 = mesh.V[tri[1] * 2 + 1];
            double x3 = mesh.V[tri[2] * 2 + 0], y3 = mesh.V[tri[2] * 2 + 1];
            double signed_twice_area = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            if (signed_twice_area < 0.0) {
                int tmp = tri[1];
                tri[1] = tri[2];
                tri[2] = tmp;
            }
        }
    }

    int NPG;
    {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        if (!(ss >> NPG)) return false;
    }

    std::vector<std::map<int, int>> PG(NPG);
    for (int i = 0; i < NPG; ++i) {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        int Npgn;
        std::string type;
        if (!(ss >> Npgn >> type)) return false;
        for (int n = 0; n < Npgn; ++n) {
            std::string line2;
            if (!std::getline(f, line2)) return false;
            std::istringstream ss2(line2);
            int n1, n2;
            if (!(ss2 >> n1 >> n2)) return false;
            n1--; n2--;
            PG[i][n1] = n2;
            PG[i][n2] = n1;
        }
    }
    f.close();

    mesh.Area.resize(mesh.Ne);
    std::vector<int> I2E;
    std::vector<double> In, In_len;
    std::vector<int> B2E;
    std::vector<double> Bn, Bn_len;
    std::vector<int> BedgeNodeOffset;
    std::vector<int> BedgeNodes;
    BedgeNodeOffset.push_back(0);
    std::map<EdgeKey, VisitedEntry> visited;

    for (int i = 0; i < mesh.Ne; ++i) {
        int* tri = &mesh.E[i * 3];
        double x1 = mesh.V[tri[0] * 2 + 0], y1 = mesh.V[tri[0] * 2 + 1];
        double x2 = mesh.V[tri[1] * 2 + 0], y2 = mesh.V[tri[1] * 2 + 1];
        double x3 = mesh.V[tri[2] * 2 + 0], y3 = mesh.V[tri[2] * 2 + 1];
        mesh.Area[i] = 0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
        if (mesh.Area[i] < 0) mesh.Area[i] = -mesh.Area[i];

        for (int j = 0; j < 3; ++j) {
            int na = tri[j], nb = tri[(j + 1) % 3];
            int end0 = na < nb ? na : nb;
            int end1 = na < nb ? nb : na;
            EdgeKey faceID(end0, end1);

            auto itB = B.find(faceID);
            if (itB != B.end()) {
                bool periodic = false;
                int pg_other0 = -1, pg_other1 = -1;
                for (size_t pg = 0; pg < PG.size(); ++pg) {
                    auto m = PG[pg];
                    if (m.count(end0) && m.count(end1)) {
                        periodic = true;
                        pg_other0 = m.at(end0);
                        pg_other1 = m.at(end1);
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
                    const BoundaryEdgeInfo& bedge = itB->second;
                    if (bedge.nodes.size() < 2) return false;
                    int bgroup = bedge.bgroup;
                    B2E.push_back(i);
                    B2E.push_back(j);
                    B2E.push_back(bgroup);
                    double nx, ny, len;
                    face_outward_normal(mesh, i, j, nx, ny, len);
                    Bn.push_back(nx);
                    Bn.push_back(ny);
                    Bn_len.push_back(len);

                    bool same_dir = (bedge.nodes[0] == na && bedge.nodes[1] == nb);
                    bool opp_dir = (bedge.nodes[0] == nb && bedge.nodes[1] == na);
                    BedgeNodes.push_back(na);
                    BedgeNodes.push_back(nb);
                    if (same_dir) {
                        for (size_t idx = 2; idx < bedge.nodes.size(); ++idx)
                            BedgeNodes.push_back(bedge.nodes[idx]);
                    } else if (opp_dir) {
                        for (int idx = (int)bedge.nodes.size() - 1; idx >= 2; --idx)
                            BedgeNodes.push_back(bedge.nodes[(size_t)idx]);
                    } else {
                        for (size_t idx = 2; idx < bedge.nodes.size(); ++idx)
                            BedgeNodes.push_back(bedge.nodes[idx]);
                    }
                    BedgeNodeOffset.push_back((int)BedgeNodes.size());
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
    mesh.BedgeNodeOffset = std::move(BedgeNodeOffset);
    mesh.BedgeNodes = std::move(BedgeNodes);

    return true;
}
