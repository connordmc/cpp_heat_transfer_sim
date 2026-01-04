import math
from math import sqrt

def Cu_thermal_con(T):
    '''Takes temperature in kelvin and returns the thermal conductivity
    using a formula from NIST for OFHC Copper'''
    a = 1.8743
    b = -0.41538
    c = -0.6018
    d = 0.13294
    e = 0.26426
    f = -0.0219
    g = -0.051276
    h = 0.0014871
    i = 0.003723
    log_10_k = (a + c* sqrt(T) + e*T + g*T*sqrt(T) + i*T*T)/(1+b*sqrt(T) + d*T + f*T*sqrt(T) + h*T*T)
    if T <= 4:
        return 10**((a + c* sqrt(4) + e*4 + g*4*sqrt(4) + i*16)/(1+b*sqrt(4) + d*4 + f*T*sqrt(4) + h*16))
    return 10**log_10_k

def Cu_specific_heat(T):
    '''
    Takes temperature in Kelvin and returns the 
    Specific Heat (J/kg-K) for OFHC Copper.
    Uses NIST polynomial for T >= 4K and a linear-cubic fit for T < 4K.
    '''
    if T < 4:
        gamma = 0.0108 
        alpha = 0.00076
        return gamma * T + alpha * (T**3)

    a, b, c = -1.91844, -0.15973, 8.61013
    d, e, f = -18.996, 21.9661, -12.7328
    g, h, i = 3.54322, -0.3797, 0 

    logT = math.log10(T)
    log10_Cp = (a + b*logT + c*(logT**2) + d*(logT**3) + 
                e*(logT**4) + f*(logT**5) + g*(logT**6) + h*(logT**7))
    
    return 10**log10_Cp
import math

import numpy as np
from scipy.integrate import quad

def Cu_resistivity(T):
    """
    Copper resistivity using Bloch-Grüneisen model.
    Returns resistivity in ohm-m.
    Valid from ~1 K to 300 K.
    """
    RRR = 50

    theta_D = 343.0  # Debye temperature for Cu (K)
    rho_300 = 15.53e-9
    rho0 = rho_300 / RRR

    def bg_integral(x):
        return (x**5) / ((np.exp(x) - 1) * (1 - np.exp(-x)))

    # Temperature-dependent phonon resistivity
    def rho_ph(T):
        upper = theta_D / T
        integral, _ = quad(bg_integral, 0, upper)
        return rho_300 * (T / theta_D)**5 * integral / (
            (273 / theta_D)**5 *
            quad(bg_integral, 0, theta_D / 273)[0]
        )

    return rho0 + rho_ph(T)

def nbti_specific_heat(T):
    # Supposes a density of 6000 kg/m^3
    if T<20: # this stopped at tc
        return 0.165714286*T+.0029*T**3
    if T>20 and T<50:
        return 7.392857143-1.401089286*T+0.098876786*T**2+0.002139286*T**3 - 0.000038929*T**4
    if 50<T and T<175:
        return -273.2142857 +14.82535714*T - 0.127910714*(T**2) + 0.000531429*(T**3) - 8.607142857e-7*(T**4)
    if T>175: # this stopped at 500K
        return 221.4285714 + 2.4475*T - 0.009225*(T**2) + 1.66e-5*(T**3) - 1.123214286e-8*(T**4)
    
def nbti_thermal_conductivity(T):
    if 10>=T: # this stopped at 1K
        return -0.000890853506962360*(T**3) + 0.016706386304553200*(T**2) - 0.044789876496699500*T + 0.068105653491378900
    if T>10 and 100>=T:
        return 2.707071295071070e-14*(T**6) - 1.089542905638570e-10*(T**5) + 7.261426646553600e-08*(T**4) - 1.768885700474560e-05 * (T**3) + 1.523577906200000e-03 * (T**2) + 1.965743220116850e-02*T + 6.411246994511720e-04
    if T>100: # data stopped at 300K
        return  -3.611887887990550e-08*T**3 + 8.451698778158840e-05*(T**2) -  1.476621990221600e-02*T + 6.425342077064450
    
def nbti_electrical_resistivity(T):
    # https://cds.cern.ch/record/527184/files/rpph096.pdf
    Tc = 9.2
    if T>Tc:
        return 5.6e-7 #ohm meters
    if T<=Tc:
        return 0



if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    T = np.linspace(0.5, 300, 2000)

    # --- Electrical Resistivity ---
    rho_cu = np.array([Cu_resistivity(t) for t in T])
    rho_nbti = np.array([nbti_electrical_resistivity(t) for t in T])

    plt.figure()
    plt.semilogy(T, rho_cu, label="OFHC Copper (RRR=50)")
    plt.semilogy(T, rho_nbti, label="NbTi")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Electrical Resistivity [Ω·m]")
    plt.title("Electrical Resistivity vs Temperature")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.show()

    # --- Thermal Conductivity ---
    k_cu = np.array([Cu_thermal_con(t) for t in T])
    k_nbti = np.array([nbti_thermal_conductivity(t) for t in T])

    plt.figure()
    plt.semilogy(T, k_cu, label="OFHC Copper")
    plt.semilogy(T, k_nbti, label="NbTi")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Thermal Conductivity [W/m·K]")
    plt.title("Thermal Conductivity vs Temperature")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.show()

    # --- Specific Heat ---
    cp_cu = np.array([Cu_specific_heat(t) for t in T])
    cp_nbti = np.array([nbti_specific_heat(t) for t in T])

    plt.figure()
    plt.loglog(T, cp_cu, label="OFHC Copper")
    plt.loglog(T, cp_nbti, label="NbTi")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific Heat [J/kg·K]")
    plt.title("Specific Heat vs Temperature")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.show()
    
    print(8960*cp_cu.min()/(k_cu.max()*2))
    print(6000*cp_nbti.min()/(k_nbti.max()*2))