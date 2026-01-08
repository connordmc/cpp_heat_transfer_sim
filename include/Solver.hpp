#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "Node.hpp"
#include "Materials.hpp"
#include <vector>

class InitialSolver{
    public:
    InitialSolver(std::vector<Node>& node, double delx);
    ~InitialSolver() = default;
    virtual double step_dynamic(double dt, double current) = 0;
    virtual int step(double dt, double current) = 0;
    std::vector<double> get_temps() const;
    const std::vector<Node>& get_nodes() const;

    protected:
    double dx;
    std::vector<Node> nodes;
    void update_material_properties();
    double qgen(double elec_res, double current);
    double qgen_base(double resistance, double current);
};

class ForwardEulerSolver : public InitialSolver{    
    public:
    ForwardEulerSolver(std::vector<Node>& node, double delx);
    double step_dynamic(double dt, double current) override;
    int step(double dt, double current) override;
};

class CrankNicolsonSolver : public InitialSolver{
    public:
    CrankNicolsonSolver(std::vector<Node>& node, double delx);
    int step(double dt, double current) override;
    double step_dynamic(double dt, double current) override; //Temporary, will implement later.

};

#endif