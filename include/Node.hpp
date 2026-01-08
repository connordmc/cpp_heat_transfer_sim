#ifndef NODE_HPP
#define NODE_HPP
#include <cstddef>
#include "Materials.hpp"

struct Node{
    Node(double T, const Material& mat, bool fixed);
    Node(double T, const Material& mat, bool fixed, double power);
    double T;

    const Material& mat;
    double k, rho_e, c, rho;
    bool fixed;

    size_t idx_k, idx_rho_e, idx_c;
    double max_power, fixed_temp; // only used on fixed nodes

    void add_max_power(double power);
};
#endif