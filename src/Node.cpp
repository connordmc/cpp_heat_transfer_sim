#include "../include/Node.hpp"

Node::Node(double Temp, const Material& m, bool fixed)
    : T(Temp),
      mat(m),
      idx_k(0),
      idx_rho_e(0),
      idx_c(0),
      rho(m.rho),
      fixed(fixed),
      max_power(0.0)
{}

Node::Node(double Temp, const Material& m, bool fixed, double power, double goal_temp)
    : T(Temp),
      mat(m),
      idx_k(0),
      idx_rho_e(0),
      idx_c(0),
      rho(m.rho),
      fixed(fixed),
      max_power(power),
      fixed_temp(goal_temp)
{}

void Node::add_max_power(double power){
  if (!fixed){
    throw std::logic_error("max_power set on non-fixed node");
  }
  max_power = power;
}

Plate::Plate(double x_center, double thickness, double initial_temp, double max_power):
  x_center(x_center), thickness(thickness), initial_temp(initial_temp), max_power(max_power){}

