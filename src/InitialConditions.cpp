#include "../include/InitialConditions.hpp"


std::vector<double> InitialConditions::InitTempsLinear(
    const std::vector<double>& plate_temps,
    const std::vector<int>& plate_indices,
    int node_count, const double& dx)
{
    std::vector<double> temps(node_count, 0);

    int curr_plate = 0;

    double x0 = plate_indices[curr_plate] * dx;
    double x1 = plate_indices[curr_plate + 1] * dx;

    double y0 = plate_temps[curr_plate];
    double y1 = plate_temps[curr_plate + 1];

    for (int i = 0; i < node_count; ++i)
    {
        if (curr_plate + 1 < (int)plate_indices.size() &&
            i > plate_indices[curr_plate + 1])
        {
            ++curr_plate;
            x0 = plate_indices[curr_plate] * dx;
            x1 = plate_indices[curr_plate + 1] * dx;

            y0 = plate_temps[curr_plate];
            y1 = plate_temps[curr_plate + 1];
        }

        if (i == plate_indices[curr_plate])
        {
            temps[i] = plate_temps[curr_plate];
            continue;
        }

        if (curr_plate + 1 >= (int)plate_indices.size())
        {
            temps[i] = plate_temps[curr_plate];
            continue;
        }

        double t = (i * dx - x0) / (x1 - x0);
        temps[i] = y0 + t * (y1 - y0);
    }

    return temps;
}
