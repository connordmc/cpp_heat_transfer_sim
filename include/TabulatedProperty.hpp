#ifndef TABULATED_PROPERTY_HPP
#define TABULATED_PROPERTY_HPP
#include <vector>


struct TabulatedProperty{
    std::vector<double> T;
    std::vector<double> vals;
    double eval(double temp, size_t& idx) const;
};

#endif