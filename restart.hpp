#ifndef RESTART_HPP
#define RESTART_HPP

#include "readgri.hpp"
#include <string>
#include <vector>

struct RestartData {
    int order = -1;
    int Ne = 0;
    int Np = 0;
    int nvars = 4;
    std::string mesh_file;
    std::vector<double> U; // layout: U[var * Ne * Np + elem * Np + basis]
};

/**
 * Write DG coefficients to a restart .dat file.
 */
bool write_restart_dat(const std::string& filepath,
                       const GriMesh& mesh,
                       int order,
                       const double* U,
                       const std::string& mesh_file);

/**
 * Read DG coefficients from a restart .dat file.
 */
bool read_restart_dat(const std::string& filepath, RestartData& data);

/**
 * Initialize target-order solution from source restart data on the same mesh.
 * Supports p-up, p-down, and same-p via nodal interpolation.
 */
bool initialize_from_restart(const RestartData& src,
                             int target_order,
                             int target_Ne,
                             double* U_target,
                             std::string& err_msg);

#endif
