#ifndef READGRI_HPP
#define READGRI_HPP

#include <vector>
#include <string>

struct GriMesh {
    int Nn;
    int Ne;
    std::vector<double> V;
    std::vector<int> E;
    std::vector<double> Area;

    int num_interior_faces;
    std::vector<int> I2E;
    std::vector<double> In;
    std::vector<double> In_len;

    int num_boundary_faces;
    std::vector<int> B2E;
    std::vector<double> Bn;
    std::vector<double> Bn_len;
    std::vector<int> BedgeNodeOffset;
    std::vector<int> BedgeNodes;

    std::vector<std::string> Bname;
};

bool read_gri(const char* fname, GriMesh& mesh);

#endif
