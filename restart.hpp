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
    std::vector<double> U;
};

bool write_restart_dat(const std::string& filepath,
                       const GriMesh& mesh,
                       int order,
                       const double* U,
                       const std::string& mesh_file);

bool read_restart_dat(const std::string& filepath, RestartData& data);

bool initialize_from_restart(const RestartData& src,
                             int target_order,
                             int target_Ne,
                             double* U_target,
                             std::string& err_msg);

bool initialize_from_restart_cross_mesh(const RestartData& src,
                                        const GriMesh& src_mesh,
                                        int target_order,
                                        const GriMesh& target_mesh,
                                        double* U_target,
                                        std::string& err_msg);

#endif
