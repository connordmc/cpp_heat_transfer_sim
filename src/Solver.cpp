#include "../include/Solver.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>

InitialSolver::InitialSolver(std::vector<Node>& node, double delx)
    :nodes(node), dx(delx)
{};

double InitialSolver::qgen(double elec_res, double current){
    return elec_res * current*current; // suppose cross sectional area = unity
}

double InitialSolver::qgen_base(double resistance, double current){
    return resistance*current*current;
}
std::vector<double> InitialSolver::get_temps() const {
    std::vector<double> temps;
    temps.reserve(nodes.size());
    for (const auto& n : nodes) {
        temps.push_back(n.T);
    }
    return temps;
}

const std::vector<Node>& InitialSolver::get_nodes() const{
    return nodes;
}

void InitialSolver::update_material_properties(){
    for (Node& n: nodes){
        n.k = n.mat.k_lookup.eval(n.T, n.idx_k);
        n.c = n.mat.cp_lookup.eval(n.T, n.idx_c);
        n.rho_e = n.mat.rho_elec_lookup.eval(n.T, n.idx_rho_e);
    }
}

ForwardEulerSolver::ForwardEulerSolver(std::vector<Node>& node, double delx) : InitialSolver(node, delx){}

int ForwardEulerSolver::step(double dt, double current){
    update_material_properties();

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
    return 0;
}


double ForwardEulerSolver::step_dynamic(double dt, double current){
    double dt_min = 1;
    update_material_properties();

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

CrankNicolsonSolver::CrankNicolsonSolver(std::vector<Node>& node, double delx) : InitialSolver(node, delx){}

int CrankNicolsonSolver::step(double dt, double current){
    double cross_sectional_area = 0.0005*0.0005*3.14159; // 1 mm diameter
    update_material_properties();
    int num_nodes = nodes.size();
    Eigen::SparseMatrix<double> A(num_nodes, num_nodes);
    Eigen::SparseMatrix<double> B(num_nodes, num_nodes);
    std::vector<Eigen::Triplet<double>> tripletsA;
    std::vector<Eigen::Triplet<double>> tripletsB;
    Eigen::VectorXd T_vec = Eigen::VectorXd::Zero(num_nodes);
    Eigen::VectorXd G = Eigen::VectorXd::Zero(num_nodes);
    Eigen::VectorXd Jacob = Eigen::VectorXd::Zero(num_nodes);
    double dTemp = 0.005;
    double Gp, Gn;
    double Jmax = 10/dt;
    double resistor_length = 0.003;

    for (int i=0; i<num_nodes; i++){
        Node& n = nodes[i];
        T_vec(i) = n.T;
        MaterialNextState next_vals = compute_gamma_at_node(n, dTemp);
        double factor = dt/(n.rho*n.c*cross_sectional_area*dx);
        double next_factor = dt/(n.rho*next_vals.next_cp*cross_sectional_area*dx);
        double prev_factor = dt/(n.rho*next_vals.prev_cp*cross_sectional_area*dx);
        double R = dt*n.k / (2*n.rho*n.c*dx*dx);
        if (i == num_nodes-1){
            tripletsA.emplace_back(Eigen::Triplet<double>(i, i, 1.0));
            //tripletsB.emplace_back(Eigen::Triplet<double>(i, i, 1.0));
            G(i) = 300.0;
            continue;
        }
        
        else if (i==0){
            tripletsA.emplace_back(Eigen::Triplet<double>(i, i, 1 + 2*R));
            tripletsB.emplace_back(Eigen::Triplet<double>(i, i, 1.0));
            G(i) = factor*qgen_base(1, current);
            Jacob(i) = 0.0;
        }
        else{

            if (n.fixed) {
                double temp_differential = n.T - n.fixed_temp;
                double smoothing = dx *100* (std::tanh(temp_differential/10));
                G(i) = factor*(qgen(n.rho_e, current) - smoothing*n.max_power);
                Gp = prev_factor*(qgen(next_vals.prev_rho_e, current) - smoothing*n.max_power);
                Gn = next_factor*(qgen(next_vals.next_rho_e, current) - smoothing*n.max_power);
            }
            else{
                G(i) = factor*qgen(n.rho_e, current);
                Gp = prev_factor*(qgen(next_vals.prev_rho_e, current));
                Gn = next_factor*(qgen(next_vals.next_rho_e, current));
            }
            if (i == num_nodes -2){
                G(i) += 600*R;
            }
            Jacob(i) = std::clamp((Gn - Gp) / (2*dTemp), -Jmax, Jmax);

            tripletsA.emplace_back(Eigen::Triplet<double>(i, i, 1 + 2*R + Jacob(i)));
            tripletsB.emplace_back(Eigen::Triplet<double>(i, i, 1 - 2*R));
            if (i > 0) {
                tripletsA.emplace_back(Eigen::Triplet<double>(i, i-1, -R));
                tripletsB.emplace_back(Eigen::Triplet<double>(i, i-1, R));
            }
            if (i < num_nodes-2) {
                tripletsA.emplace_back(Eigen::Triplet<double>(i, i+1, -R));
                tripletsB.emplace_back(Eigen::Triplet<double>(i, i+1, R));
            }
        }

    }
    A.setFromTriplets(tripletsA.begin(), tripletsA.end());
    B.setFromTriplets(tripletsB.begin(), tripletsB.end());
    Eigen::VectorXd b_RHS = B*T_vec + G + Jacob.cwiseProduct(T_vec);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed" << std::endl;
        return 1;
    }

    Eigen::VectorXd T_new = solver.solve(b_RHS);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed" << std::endl;
        return 1;
    }

    for (int i=0; i<num_nodes-1; i++){
        Node& n = nodes[i];
        n.T = T_new[i];
    }
    nodes[num_nodes-1].T = 300.0; // The outside, infinite drain
    return 0;
}

double CrankNicolsonSolver::step_dynamic(double dt, double current){
    throw std::runtime_error("Error: Not yet implemented");
}

MaterialNextState CrankNicolsonSolver::compute_gamma_at_node(Node& node, double dT){
    /*
    Returns a vector of the form <rho_e +1, c_p +1, rho_e -1, c_p -1>
    */
    const Material& mat = node.mat;
    MaterialNextState next_state;
    double T = node.T;
    next_state.next_rho_e = mat.rho_elec_lookup.eval(T + dT, node.idx_rho_e);
    next_state.next_cp = mat.cp_lookup.eval(T + dT, node.idx_c);
    next_state.prev_rho_e = mat.rho_elec_lookup.eval(T - dT, node.idx_rho_e);
    next_state.prev_cp = mat.cp_lookup.eval(T - dT, node.idx_c);
    return next_state;
}