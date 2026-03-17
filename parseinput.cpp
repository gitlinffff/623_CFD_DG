#include "parseinput.hpp"
#include <fstream>
#include <sstream>

bool read_input_file(const std::string& file, InputParams& p)
{
    std::ifstream f(file);
    if (!f) return false;

    std::string line, key, val;

    while (std::getline(f, line)) {
        if (line.empty() || line[0]=='#') continue;

        std::istringstream iss(line);
        if (std::getline(iss,key,'=') && std::getline(iss,val)) {

            if (key=="order")               p.order = std::stoi(val);
            else if (key=="CFL")            p.CFL = std::stod(val);
            else if (key=="gri_file")       p.gri_file = val;
            else if (key=="max_iter")       p.max_iter = std::stoi(val);
            else if (key=="print_interval") p.print_interval = std::stoi(val);
            else if (key=="use_local_dt")   p.use_local_dt = std::stoi(val);
            else if (key=="inflow_perturb") p.in_ptb = std::stoi(val);
						else if (key=="flux")           p.flux = val;
        }
    }

    return true;
}
