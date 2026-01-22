#include "../include/Tests.hpp"
#include "../include/Solver.hpp"
#include "../include/Node.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace Tests{
    void check_initial_conditions(const std::vector<double>& temps, int step){
        for (size_t i = 0; i < temps.size(); i++) {
            if (!std::isfinite(temps[i]))
                throw std::runtime_error("NaN/Inf temperature at node " + std::to_string(i) + " on step " + std::to_string(step));

            if (temps[i] < 0.0){
                throw std::runtime_error("Negative temperature at node " + std::to_string(i) + " on step " + std::to_string(step));
            }
            if (temps[i] > 10000){
                throw std::runtime_error("Temperature exceeds physical expectations at node " + std::to_string(i)
                + " on step " + std::to_string(step));
            }
        }
    }
    void check_for_tc(const std::vector<double>& temps, int step){
        for (size_t i = 0; i<temps.size(); i++){
            if (temps[i]<9.2){
                std::cout << "Temperature below Tc occurs at step " << std::to_string(i) << " on step " << std::to_string(step) << std::endl;
            }
        }
    }

    void check_maximum_power(const Node& n, double area, double dx, double dt, double safety_ratio){
        if (n.max_power <= 0.0){
            return;
        }

        if (n.T <= 0.0) {
            throw std::runtime_error("Invalid node temperature (T <= 0)");
        }

        const double E_node = n.rho * n.c * area * dx * n.T;
        const double safe_power = safety_ratio * E_node / dt;
        const double safe_dt = safety_ratio * E_node / n.max_power;
        if (n.max_power > safe_power){
            throw std::runtime_error("Power of " + std::to_string(n.max_power) + " exceeds maximum, lower to " + std::to_string(safe_power)
            + " or reduce dt to " + std::to_string(safe_dt) + "E_node = " + std::to_string(E_node));
        }

    }
}