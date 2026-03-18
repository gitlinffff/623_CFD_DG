#ifndef WRITE_VTU_HPP
#define WRITE_VTU_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include <string>

void write_solution_vtu(const GriMesh& mesh, const double* U, int order,
                        const ProblemParams& params, const std::string& filepath);

void write_pvd_manifest(const std::string& case_dir, const std::vector<double>& times);
#endif
