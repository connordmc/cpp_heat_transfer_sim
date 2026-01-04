#include "../include/Node.hpp"

Node::Node(double Temp, const Material& m, bool fixed)
    : T(Temp),
      mat(m),
      idx_k(0),
      idx_rho_e(0),
      idx_c(0),
      rho(m.rho),
      fixed(fixed)
{}

