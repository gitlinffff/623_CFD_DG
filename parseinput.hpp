#pragma once
#include <string>

struct InputParams {
    int order = 0;
    double CFL = 0.5;
    std::string gri_file = "mesh/global_refine_0.gri";
    int max_iter = 10000;
    int print_interval = 100;
    bool use_local_dt = false;
		std::string flux = "HLLE";
    std::string restart_file = "";
};

bool read_input_file(const std::string& file, InputParams& p);
