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

Node::Node(double Temp, const Material& m, bool fixed, double power)
    : T(Temp),
      mat(m),
      idx_k(0),
      idx_rho_e(0),
      idx_c(0),
      rho(m.rho),
      fixed(fixed),
      max_power(power)
{}

void Node::add_max_power(double power){
  if (!fixed){
    throw std::logic_error("max_power set on non-fixed node");
  }
  max_power = power;
}

