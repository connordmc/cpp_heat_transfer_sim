#include <iostream>
#include <vector>
#include <array>
#include <cstddef>
#include "../include/Solver.hpp"
#include "../include/MaterialLoad.hpp"
#include "H5Cpp.h"
#include <chrono>
#include <iomanip>

int main(){
    // Setup
    constexpr std::size_t NumPlates = 7;
    constexpr std::size_t NumNodes  = 1000;

    Material cu   = MaterialLoad::load_copper("data");
    Material nbti = MaterialLoad::load_nbti("data");

    constexpr double nbti_min = 0.4;
    constexpr double nbti_max = 0.6;

    const std::array<double, NumPlates> relative_pos_plates =
        {0.0, 0.15, 0.15, 0.15, 0.15, 0.40, 0.40};

    const std::array<double, NumPlates> initial_plate_temps =
        {0.0, 0.01, 0.1, 1.0, 4.0, 50.0, 300.0};

    std::array<double, NumPlates> pos_plates{};
    std::array<int,   NumPlates> plate_indices{};

    double length = 0.0f;
    for (std::size_t i = 0; i < NumPlates; ++i) {
        length += relative_pos_plates[i];
        pos_plates[i] = length;
    }

    const double dx = length / static_cast<double>(NumNodes);

    for (std::size_t i = 0; i < NumPlates; ++i) {
        plate_indices[i] =
            static_cast<int>(std::round(pos_plates[i] / dx));
    }

    std::vector<Node> nodes;
    nodes.reserve(NumNodes);

    std::vector<char> is_fixed_node(NumNodes, 0);
    std::vector<double> initial_temps(NumNodes, 300.0f);

    for (std::size_t p = 0; p < NumPlates; ++p) {
        if (initial_plate_temps[p] != 0.0f) {
            int idx = plate_indices[p];
            if (idx >= 0 && idx < static_cast<int>(NumNodes)) {
                is_fixed_node[idx] = 1;
                initial_temps[idx] = initial_plate_temps[p];
            }
        }
    }

    for (std::size_t i = 0; i < NumNodes; ++i) {
        double x = i * dx;

        const Material& mat =
            (x >= nbti_min && x <= nbti_max) ? nbti : cu;

        nodes.emplace_back(initial_temps[i], mat, is_fixed_node[i]!=0);
    }
    Solver heat_sim(nodes, dx);
    double dt = dx*dx*0.0167*0.100;
    constexpr std::size_t total_time = 900; // seconds
    constexpr double current = 0.0; // amps
    constexpr int save_interval = 1000;
    int total_steps = static_cast<int>(std::round(total_time/dt));
    int num_snapshots = total_steps / save_interval + 2;

    H5::H5File file("simulation_output.h5", H5F_ACC_TRUNC);
    hsize_t dims[2] = {static_cast<hsize_t>(num_snapshots), static_cast<hsize_t>(NumNodes)};
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset = file.createDataSet("temperatures", H5::PredType::NATIVE_DOUBLE, dataspace);
    {
        H5::Attribute attr_dx = dataset.createAttribute("dx", H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
        attr_dx.write(H5::PredType::NATIVE_DOUBLE, &dx);

        H5::Attribute attr_dt = dataset.createAttribute("dt", H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
        attr_dt.write(H5::PredType::NATIVE_DOUBLE, &dt);

        hsize_t dims_plate[1] = {plate_indices.size()};
        H5::DataSpace plate_space(1, dims_plate);
        H5::Attribute attr_plate = dataset.createAttribute("plate_indices",
            H5::PredType::NATIVE_INT,
            plate_space);
        attr_plate.write(H5::PredType::NATIVE_INT, plate_indices.data());
    }

    // Simulation loop
    int snapshot_index = 0;
    std::vector<double> temp_row(NumNodes);
    auto start_time = std::chrono::steady_clock::now();

    for (int step = 0; step <= total_steps; ++step) {
        dt = heat_sim.step(dt, current);
        if (step % (total_steps/10000) == 0){
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();
            double remaining = elapsed * (total_steps - step) / step;
            double percent = 100.00f * step/total_steps;
            std::cout << "\r" << std::fixed << std::setprecision(2) << percent << "% done, ETA: " << remaining << " s" << std::flush;
        }

        if (step % save_interval == 0) {
            temp_row = heat_sim.get_temps();

            hsize_t start[2] = {static_cast<hsize_t>(snapshot_index), 0};
            hsize_t count[2] = {1, static_cast<hsize_t>(NumNodes)};
            H5::DataSpace file_space = dataset.getSpace();
            file_space.selectHyperslab(H5S_SELECT_SET, count, start);
            hsize_t mem_dims[2] = {1, static_cast<hsize_t>(NumNodes)};
            H5::DataSpace mem_space(2, mem_dims);
            dataset.write(temp_row.data(), H5::PredType::NATIVE_DOUBLE, mem_space, file_space);
            ++snapshot_index;
        }
    }
    temp_row = heat_sim.get_temps();

    hsize_t start[2] = {static_cast<hsize_t>(snapshot_index), 0};
    hsize_t count[2] = {1, static_cast<hsize_t>(NumNodes)};

    H5::DataSpace file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    hsize_t mem_dims[2] = {1, static_cast<hsize_t>(NumNodes)};
    H5::DataSpace mem_space(2, mem_dims);
    dataset.write(temp_row.data(), H5::PredType::NATIVE_DOUBLE, mem_space, file_space);

    return 0;
}
