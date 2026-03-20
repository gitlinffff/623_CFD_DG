#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "readgri.hpp"
#include <Eigen/Dense>
#include <vector>

struct GeomEval {
    double x;
    double y;
    Eigen::Matrix2d J;
    double detJ;
};

int geom_num_nodes(int q);
int elem_geom_order(const GriMesh& mesh, int elem);

void faceRefCoords(int face, double t, double& xi, double& eta);
void faceRefDerivatives(int face, double& dxi_dt, double& deta_dt);

void eval_geometry(const GriMesh& mesh, int elem, double xi, double eta, GeomEval& g);
bool phys_to_ref(const GriMesh& mesh, int elem, double x, double y, double& xi, double& eta);

double face_metric_normal_point(const GriMesh& mesh, int elem, int face, double t,
                                double n_out[2], double& x, double& y);

const std::vector<double>& get_geometric_element_areas(const GriMesh& mesh);
double get_geometric_element_area(const GriMesh& mesh, int elem);

#endif
