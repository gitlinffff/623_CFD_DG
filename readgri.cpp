#include "readgri.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <cmath>
#include <cstring>

namespace {

using EdgeKey = std::pair<int, int>;

struct VisitedEntry {
    int elemL;   /* 0-based */
    int faceL;   /* 0-based local face */
    bool aligned;
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

}  // namespace

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

    std::map<EdgeKey, std::string> B;
    mesh.Bname.clear();
    for (int i = 0; i < NB; ++i) {
        std::string line;
        if (!std::getline(f, line)) return false;
        std::istringstream ss(line);
        int Nb;
        std::string name;
        if (!(ss >> Nb)) return false;
        ss >> name;  /* skip one token, then get name - format is "Nb ? name" */
        if (!(ss >> name)) return false;
        mesh.Bname.push_back(name);
        for (int n = 0; n < Nb; ++n) {
            std::string line2;
            if (!std::getline(f, line2)) return false;
            std::istringstream ss2(line2);
            int n1, n2;
            if (!(ss2 >> n1 >> n2)) return false;
            n1--; n2--;
            int a = n1 < n2 ? n1 : n2;
            int b = n1 < n2 ? n2 : n1;
            B[EdgeKey(a, b)] = name;
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
                        double ex = mesh.V[end1 * 2 + 0] - mesh.V[end0 * 2 + 0];
                        double ey = mesh.V[end1 * 2 + 1] - mesh.V[end0 * 2 + 1];
                        double len = norm2(ex, ey);
                        if (len < 1e-20) len = 1e-20;
                        ex /= len; ey /= len;
                        double nx = ey, ny = -ex;
                        double ox = mesh.V[o1 * 2 + 0] - mesh.V[o0 * 2 + 0];
                        double oy = mesh.V[o1 * 2 + 1] - mesh.V[o0 * 2 + 1];
                        bool is_antiparallel = (ex * ox + ey * oy) < 0;
                        bool flip = !ve.aligned;
                        if (is_antiparallel) flip = !flip;
                        if (flip) { nx = -nx; ny = -ny; }
                        In.push_back(nx);
                        In.push_back(ny);
                        In_len.push_back(len);
                        visited.erase(itV);
                    } else {
                        visited[faceID] = { i, j, na < nb };
                    }
                } else {
                    int bgroup = 1;
                    for (size_t k = 0; k < mesh.Bname.size(); ++k) {
                        if (mesh.Bname[k] == itB->second) { bgroup = (int)k + 1; break; }
                    }
                    B2E.push_back(i);
                    B2E.push_back(j);
                    B2E.push_back(bgroup);
                    double ex = mesh.V[nb * 2 + 0] - mesh.V[na * 2 + 0];
                    double ey = mesh.V[nb * 2 + 1] - mesh.V[na * 2 + 1];
                    double len = norm2(ex, ey);
                    if (len < 1e-20) len = 1e-20;
                    ex /= len; ey /= len;
                    Bn.push_back(ey);
                    Bn.push_back(-ex);
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
                    double ex = mesh.V[end1 * 2 + 0] - mesh.V[end0 * 2 + 0];
                    double ey = mesh.V[end1 * 2 + 1] - mesh.V[end0 * 2 + 1];
                    double len = norm2(ex, ey);
                    if (len < 1e-20) len = 1e-20;
                    ex /= len; ey /= len;
                    double nx = ey, ny = -ex;
                    if (!ve.aligned) { nx = -nx; ny = -ny; }
                    In.push_back(nx);
                    In.push_back(ny);
                    In_len.push_back(len);
                    visited.erase(itV);
                } else {
                    visited[faceID] = { i, j, na < nb };
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
