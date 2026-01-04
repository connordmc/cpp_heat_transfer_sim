#ifndef MATERIALS_HPP
#define MATERIALS_HPP
#include "TabulatedProperty.hpp"

struct Material{
    TabulatedProperty k_lookup;
    TabulatedProperty rho_elec_lookup;
    TabulatedProperty cp_lookup;
    double rho;
};
#endif