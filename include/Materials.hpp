#ifndef MATERIALS_HPP
#define MATERIALS_HPP
#include "TabulatedProperty.hpp"

struct Material{
    TabulatedProperty k_lookup;
    TabulatedProperty rho_elec_lookup;
    TabulatedProperty cp_lookup;
    double rho;
};

struct MaterialNextState{
    double rho, next_cp, prev_cp, next_rho_e, prev_rho_e;
};
#endif