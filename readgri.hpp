#ifndef READGRI_HPP
#define READGRI_HPP

#include <vector>
#include <string>

/** Mesh and face data from .gri (Gambit); all indices 0-based for C++. */
struct GriMesh {
    int Nn;
    int Ne;
    std::vector<double> V;        /* Nn*2: x,y per vertex */
    std::vector<int> E;           /* Ne*3: corner-node connectivity (triangle vertices) */
    std::vector<int> ElemGeomOrder;   /* Ne: geometric order q per element (typically 1 or 2) */
    std::vector<int> ElemGeomOffset;  /* Ne+1: offsets into ElemGeomConn */
    std::vector<int> ElemGeomConn;    /* concatenated geometric node connectivity in shapeL ordering */
    std::vector<double> Area;     /* Ne: area per element */

    int num_interior_faces;
    std::vector<int> I2E;          /* num_interior_faces * 4: [elemL, faceL, elemR, faceR] */
    std::vector<double> In;       /* num_interior_faces * 2: unit normal (L to R), use &In[2*i] as n[2] */
    std::vector<double> In_len;   /* num_interior_faces: face length */

    int num_boundary_faces;
    std::vector<int> B2E;          /* num_boundary_faces * 3: [elem, face, bgroup] */
    std::vector<double> Bn;       /* num_boundary_faces * 2: unit outward normal, use &Bn[2*i] as n[2] */
    std::vector<double> Bn_len;   /* num_boundary_faces: face length */

    std::vector<std::string> Bname; /* boundary names; bgroup in B2E is 1-based index into Bname */
};

/** Load .gri file and build face tables (normals, lengths). Returns true on success. */
bool read_gri(const char* fname, GriMesh& mesh);

#endif
