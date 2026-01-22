#ifndef TESTS_HPP
#define TESTS_HPP
#include <vector>

class Node;
class InitialSolver;

namespace Tests{
    void check_initial_conditions(const std::vector<double>& temps, int step);
    void check_for_tc(const std::vector<double>& temps, int step);
    void check_maximum_power(const Node& n, double area, double dx, double dt, double safety_ratio = 0.8);
}

#endif