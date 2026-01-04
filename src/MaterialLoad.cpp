#include "../include/MaterialLoad.hpp"

void MaterialLoad::load_csv_to_prop(TabulatedProperty& prop, const std::string& filepath){
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open material file: " + filepath);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string t_str, v_str;

        if (std::getline(ss, t_str, ',') && std::getline(ss, v_str, ',')) {
            prop.T.push_back(std::stod(t_str));
            prop.vals.push_back(std::stod(v_str));
        }
    }
}

Material MaterialLoad::load_copper(const std::string& data_dir) {
    Material cu;
    cu.rho = 8960.0; // Density kg/m^3
    load_csv_to_prop(cu.k_lookup, data_dir + "/cu_k.csv");
    load_csv_to_prop(cu.cp_lookup, data_dir + "/cu_cp.csv");
    load_csv_to_prop(cu.rho_elec_lookup, data_dir + "/cu_rho_e.csv");
    return cu;
}

Material MaterialLoad::load_nbti(const std::string& data_dir) {
    Material nbti;
    nbti.rho = 6000.0; // Density kg/m^3
    load_csv_to_prop(nbti.k_lookup, data_dir + "/nbti_k.csv");
    load_csv_to_prop(nbti.cp_lookup, data_dir + "/nbti_cp.csv");
    load_csv_to_prop(nbti.rho_elec_lookup, data_dir + "/nbti_rho_e.csv");
    return nbti;
}