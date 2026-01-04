#include "../include/TabulatedProperty.hpp"
#include <algorithm>


double TabulatedProperty::eval(double temp, size_t& idx) const{
    if (temp<T.front()) return vals.front();
    if (temp>T.back()) return vals.back();

    auto it = std::lower_bound(T.begin(), T.end(), temp);
    idx = std::distance(T.begin(), it) -1;
    double T0 = T[idx], T1 = T[idx+1];
    double v0 = vals[idx], v1 = vals[idx+1];

    double w = (v1 - v0)/(T1 - T0);
    return v0 + (temp - T0)*w; 
}