#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "Node.hpp"
#include "Materials.hpp"
#include <vector>

class Solver{
    public:
    Solver(std::vector<Node>& node, double delx);
    double step(double dt, double current);
    std::vector<double> get_temps() const;
    const std::vector<Node>& get_nodes() const;

    private:
    double dx;
    std::vector<Node> nodes;
};

double qgen(double elec_res, double current);
double qgen_base(double resistance, double current);
#endif