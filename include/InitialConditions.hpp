#ifndef INT_CON_HPP
#define INT_CON_HPP
#include <vector>

namespace InitialConditions{
    std::vector<double> InitTempsLinear(
        const std::vector<double>& plate_temps,
        const std::vector<int>& plate_indices,
        int node_count, const double& dx);
}

#endif