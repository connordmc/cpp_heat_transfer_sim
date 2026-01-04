#include "../include/Solver.hpp"

Solver::Solver(std::vector<Node>& node, double delx)
    :nodes(node), dx(delx)
{};

double Solver::step(double dt, double current){
    double dt_min = 1;
    for (Node& n: nodes){
        n.k = n.mat.k_lookup.eval(n.T, n.idx_k);
        n.c = n.mat.cp_lookup.eval(n.T, n.idx_c);
        n.rho_e = n.mat.rho_elec_lookup.eval(n.T, n.idx_rho_e);
    }

    std::vector<double> T_old(nodes.size());
    for (int i = 0; i < nodes.size(); i++){
        T_old[i] = nodes[i].T;
    }

    for (int i=1; i <nodes.size() - 1; i++ ){
        Node& n = nodes[i];

        
        if (n.fixed == false){
            const Node& n_left = nodes[i-1];
            const Node& n_right = nodes[i+1];
            double k_left = 2*n_left.k*n.k/(n_left.k+n.k) ;
            double k_right = 2*n_right.k*n.k/(n_right.k+n.k);
            double left_flux = k_left*(T_old[i-1] - T_old[i])/dx;
            double right_flux = k_right*(T_old[i+1] - T_old[i])/dx;
            double scale = dt / (n.rho * n.c);
            n.T = T_old[i] + scale*((left_flux + right_flux)/dx
            + qgen(n.rho_e, current));
            double local_safe_dt = 0.7*(n.rho * n.c*dx*dx)/(2*n.k);

            if (local_safe_dt < dt_min){
                dt_min = local_safe_dt;
            }
        }

    }
    Node& n0 = nodes[0];
    // We treat the bottom node as a simple resistor with set resistance 1 ohm.
    if (!n0.fixed) {
        double scale = dt/(n0.rho*n0.c);
        n0.T = T_old[0] + scale*(n0.k*2.0*(T_old[1]-T_old[0])/(dx*dx) + qgen_base(1.0, current));
    }
    else{
        n0.T = T_old[0];
    }
    nodes[nodes.size()-1].T = T_old[nodes.size()-1];
    return dt_min;
}

double qgen(double elec_res, double current){
    return elec_res * current*current; // suppose cross sectional area = unity
}

double qgen_base(double resistance, double current){
    return resistance*current*current;
}
std::vector<double> Solver::get_temps() const {
    std::vector<double> temps;
    temps.reserve(nodes.size());
    for (const auto& n : nodes) {
        temps.push_back(n.T);
    }
    return temps;
}

const std::vector<Node>& Solver::get_nodes() const{
    return nodes;
}