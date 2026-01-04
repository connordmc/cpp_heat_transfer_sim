import os
from generate_mats import (Cu_thermal_con, Cu_resistivity, Cu_specific_heat,
                           nbti_electrical_resistivity, nbti_specific_heat, nbti_thermal_conductivity)

def save_data(file_name, func, T_min = 4, T_max = 500, points = 1000):
    '''
    Generates T and property values then stores into a csv
    '''
    script_dir = os.path.dirname(os.path.abspath(__file__))

    base_dir = os.path.dirname(script_dir) 
    data_dir = os.path.join(base_dir, 'data')

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    filepath = os.path.join(data_dir, f"{file_name}.csv")
    step = (T_max - T_min)/points
    with open(filepath, 'w') as f:
        for i in range(points):
            t = T_min + step*i
            try:
                val = func(t)
                f.write(f"{t:.10f},{val:.10e}\n")
            except Exception as e:
                print(f"Error calculating {file_name} at T={t}: {e}")

if __name__ == "__main__":
    materials_to_export = {
        "cu_k": Cu_thermal_con,
        "cu_cp": Cu_specific_heat,
        "cu_rho_e": Cu_resistivity,
        "nbti_k": nbti_thermal_conductivity,
        "nbti_cp": nbti_specific_heat,
        "nbti_rho_e": nbti_electrical_resistivity
    }

    for name, func in materials_to_export.items():
        save_data(name, func)