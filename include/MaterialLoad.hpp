#include <fstream>
#include "Materials.hpp"
#include <sstream>

namespace MaterialLoad {
    void load_csv_to_prop(TabulatedProperty& prop, const std::string& filepath);
    Material load_copper(const std::string& data_dir);
    Material load_nbti(const std::string& data_dir);
}