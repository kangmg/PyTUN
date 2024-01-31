try:
    from pointgroup import PointGroup
    pointgroup_install = 'YES'
except ImportError:
    print("WARRNING : Module 'pointgroup' not found. But it is not nessesary.\n")
    print("pointgourp module used for auto-estimation of RXN Symmetry Number.")
    print("Your can install module for auto via 'pip install pointgroup'")
    print("")
    pointgroup_install = 'NO'

import numpy as np
# from decimal import Decimal
import sys, math, time
from array import array

# PHYSICAL CONSTANTS
GAS_CONSTANT = 8.3144621
PLANCK_CONSTANT = 6.62606957e-34
BOLTZMANN_CONSTANT = 1.3806488e-23
SPEED_OF_LIGHT = 2.99792458e10
AVOGADRO_CONSTANT = 6.0221415e23
AMU_to_KG = 1.66053886E-27
autokcal = 627.509541
kjtokcal = 4.184
atmos = 101.325
PI = 3.14159265359
k = 3.1668114E-6  # Boltzmann Constant in atomic units

# 10-point Gauss-Legendre Quadrature abscissa and weight (exact solution for up to 21st order polynomial)
x = array('d', [-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941, -0.1488743390, 0.1488743390, 0.4333953941,
                0.6794095683, 0.8650633667, 0.9739065285])
w = array('d', [0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247, 0.2955242247, 0.2692667193,
                0.2190863625, 0.1494513492, 0.0666713443])

# more symmetry number to be add
rot_symmetry_number = {'C1' : 1,
                       'Cs' : 1,
                       'C2' : 2,
                       'C2v': 2,
                       'C3v': 3,
                       'D2h': 4,
                       'D3h': 6,
                       'D5h': 10,
                       'Dinfh' : 2,
                       'D3d': 6,
                       'Td' : 12,
                       'Oh' : 24}



# Read gaussian output for the level of theory and basis set used
def level_of_theory(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    level = "none"
    bs = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().find('\\Freq\\') > -1:
            if len(inlines[i].strip().split("\\")) > 5:
                level = (inlines[i].strip().split("\\")[4])
                bs = (inlines[i].strip().split("\\")[5])
    return level + "/" + bs


# Read gaussian output for the Final SCF Energy
def final_scf_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    scf = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('SCF Done:'):
            scf = (inlines[i].strip().split()[4])
    return scf


# Read gaussian output for Zero Point Energy
def zero_point_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    zpe = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Zero-point correction'):
            zpe = (inlines[i].strip().split()[2])
    return zpe


# Read gaussian output for Thermal Correction to Gibbs Free Energy
def g_corr(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    gcorr = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Thermal correction to Gibbs Free Energy='):
            gcorr = (inlines[i].strip().split()[6])
    return gcorr


# Read gaussian output for Temperature
def temperature(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    temperature = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Temperature'):
            temperature = (inlines[i].strip().split()[1])
    return temperature


# Read gaussian output for reduced mass of normal mode of transition
def reduced_mass(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    mu = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Red. masses'):
            if mu == "none":
                mu = (inlines[i].strip().split()[3])
    return mu


# Read gaussian output for force constant of normal mode of transition
def force_constant(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    fc = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Frc consts'):
            if fc == "none":
                fc = (inlines[i].strip().split()[3])
    return fc


# Parameters B, ALPHA, a, b, d of Eckart Potential
def Bee(V_max, V_r, V_p):
    bee = (V_max ** 0.5 + (V_max - (V_p - V_r)) ** 0.5) ** 2
    return bee


def ALPHA(B, F_s, V_max, V_r, V_p):
    alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
    return alpha


def A(E, mu, ALPHA):
    a = 2 * PI * (2 * mu * E) ** 0.5 / ALPHA
    return a


def B(E, mu, V_p, V_r, ALPHA):
    b = 2 * PI * (2 * mu * (E - (V_p - V_r))) ** 0.5 / ALPHA
    return b


def D(bee, mu, ALPHA):
    d = 2 * PI * abs((2 * mu * bee - (ALPHA / 2) ** 2)) ** 0.5 / ALPHA
    return d


# Calculation of Transmission Probabilty of Eckart Potential
def T(a, b, d):
    T = (math.cosh(a + b) - math.cosh(a - b)) / (math.cosh(a + b) + math.cosh(d))
    return T


# Calculation of SINH function of Kappa
def S(V_max, E):
    S = math.sinh(((V_max - E)) / (TEMPERATURE * k))
    return S

# Read gaussian output for the level of theory and basis set used
def level_of_theory(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    level = "none"
    bs = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().find('\\Freq\\') > -1:
            if len(inlines[i].strip().split("\\")) > 5:
                level = (inlines[i].strip().split("\\")[4])
                bs = (inlines[i].strip().split("\\")[5])
    return level + "/" + bs


# Read gaussian output for the Final SCF Energy
def final_scf_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    scf = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('SCF Done:'):
            scf = (inlines[i].strip().split()[4])
    return scf


# Read gaussian output for Zero Point Energy
def zero_point_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    zpe = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Zero-point correction'):
            zpe = (inlines[i].strip().split()[2])
    return zpe


# Read gaussian output for Thermal Correction to Gibbs Free Energy
def g_corr(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    gcorr = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Thermal correction to Gibbs Free Energy='):
            gcorr = (inlines[i].strip().split()[6])
    return gcorr


# Read gaussian output for Temperature
def temperature(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    temperature = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Temperature'):
            temperature = (inlines[i].strip().split()[1])
    return temperature


# Read gaussian output for reduced mass of normal mode of transition
def reduced_mass(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    mu = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Red. masses'):
            if mu == "none":
                mu = (inlines[i].strip().split()[3])
    return mu


# Read gaussian output for force constant of normal mode of transition
def force_constant(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    fc = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Frc consts'):
            if fc == "none":
                fc = (inlines[i].strip().split()[3])
    return fc


# Parameters B, ALPHA, a, b, d of Eckart Potential
def Bee(V_max, V_r, V_p):
    bee = (V_max ** 0.5 + (V_max - (V_p - V_r)) ** 0.5) ** 2
    return bee


def ALPHA(B, F_s, V_max, V_r, V_p):
    alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
    return alpha


def A(E, mu, ALPHA):
    a = 2 * PI * (2 * mu * E) ** 0.5 / ALPHA
    return a


def B(E, mu, V_p, V_r, ALPHA):
    b = 2 * PI * (2 * mu * (E - (V_p - V_r))) ** 0.5 / ALPHA
    return b


def D(bee, mu, ALPHA):
    d = 2 * PI * abs((2 * mu * bee - (ALPHA / 2) ** 2)) ** 0.5 / ALPHA
    return d


# Calculation of Transmission Probabilty of Eckart Potential
def T(a, b, d):
    T = (math.cosh(a + b) - math.cosh(a - b)) / (math.cosh(a + b) + math.cosh(d))
    return T


# Calculation of SINH function of Kappa
def S(V_max, E):
    S = math.sinh(((V_max - E)) / (TEMPERATURE * k))
    return S

# Read gaussian output for the level of theory and basis set used
def level_of_theory_freq(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    level = "none"
    bs = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().find('\\Freq\\') > -1:
            if len(inlines[i].strip().split("\\")) > 5:
                level = (inlines[i].strip().split("\\")[4])
                bs = (inlines[i].strip().split("\\")[5])
    return level + "/" + bs

def level_of_theory_SP(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    level = "none"
    bs = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().find('\\SP\\') > -1:
            if len(inlines[i].strip().split("\\")) > 5:
                level = (inlines[i].strip().split("\\")[4])
                bs = (inlines[i].strip().split("\\")[5])
    return level + "/" + bs


# Read gaussian output for the Final SCF Energy
def final_scf_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    scf = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('SCF Done:'):
            scf = (inlines[i].strip().split()[4])
    return scf


# Read gaussian output for Zero Point Energy
def zero_point_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    zpe = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Zero-point correction'):
            zpe = (inlines[i].strip().split()[2])
    return zpe


# Read gaussian output for Thermal Correction to Gibbs Free Energy
def g_corr(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    gcorr = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Thermal correction to Gibbs Free Energy='):
            gcorr = (inlines[i].strip().split()[6])
    return gcorr


# Read gaussian output for Temperature
def temperature(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    temperature = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Temperature'):
            temperature = (inlines[i].strip().split()[1])
    return temperature


# Read gaussian output for reduced mass of normal mode of transition
def reduced_mass(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    mu = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Red. masses'):
            if mu == "none":
                mu = (inlines[i].strip().split()[3])
    return mu


# Read gaussian output for force constant of normal mode of transition
def force_constant(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    fc = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Frc consts'):
            if fc == "none":
                fc = (inlines[i].strip().split()[3])
    return fc


# Parameters B, ALPHA, a, b, d of Eckart Potential
def Bee(V_max, V_r, V_p):
    bee = (V_max ** 0.5 + (V_max - (V_p - V_r)) ** 0.5) ** 2
    return bee


def ALPHA(B, F_s, V_max, V_r, V_p):
    alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
    return alpha


def A(E, mu, ALPHA):
    a = 2 * PI * (2 * mu * E) ** 0.5 / ALPHA
    return a


def B(E, mu, V_p, V_r, ALPHA):
    b = 2 * PI * (2 * mu * (E - (V_p - V_r))) ** 0.5 / ALPHA
    return b


def D(bee, mu, ALPHA):
    d = 2 * PI * abs((2 * mu * bee - (ALPHA / 2) ** 2)) ** 0.5 / ALPHA
    return d


# Calculation of Transmission Probabilty of Eckart Potential
def T(a, b, d):
    T = (math.cosh(a + b) - math.cosh(a - b)) / (math.cosh(a + b) + math.cosh(d))
    return T


# Calculation of SINH function of Kappa
def S(V_max, E):
    S = math.sinh(((V_max - E)) / (TEMPERATURE * k))
    return S

# Read gaussian output for the level of theory and basis set used
def level_of_theory(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    level = "none"
    bs = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().find('\\Freq\\') > -1:
            if len(inlines[i].strip().split("\\")) > 5:
                level = (inlines[i].strip().split("\\")[4])
                bs = (inlines[i].strip().split("\\")[5])
    return level + "/" + bs


# Read gaussian output for the Final SCF Energy
def final_scf_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    scf = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('SCF Done:'):
            scf = (inlines[i].strip().split()[4])
    return scf


# Read gaussian output for Zero Point Energy
def zero_point_energy(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    zpe = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Zero-point correction'):
            zpe = (inlines[i].strip().split()[2])
    return zpe


# Read gaussian output for Thermal Correction to Gibbs Free Energy
def g_corr(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    gcorr = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Thermal correction to Gibbs Free Energy='):
            gcorr = (inlines[i].strip().split()[6])
    return gcorr


# Read gaussian output for Temperature
def temperature(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    temperature = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Temperature'):
            temperature = (inlines[i].strip().split()[1])
    return temperature


# Read gaussian output for reduced mass of normal mode of transition
def reduced_mass(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    mu = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Red. masses'):
            if mu == "none":
                mu = (inlines[i].strip().split()[3])
    return mu


# Read gaussian output for force constant of normal mode of transition
def force_constant(file):
    g09_output = open(file, 'r')
    inlines = g09_output.readlines()
    fc = "none"
    for i in range(0, len(inlines)):
        if inlines[i].strip().startswith('Frc consts'):
            if fc == "none":
                fc = (inlines[i].strip().split()[3])
    return fc


# Parameters B, ALPHA, a, b, d of Eckart Potential
def Bee(V_max, V_r, V_p):
    bee = (V_max ** 0.5 + (V_max - (V_p - V_r)) ** 0.5) ** 2
    return bee


def ALPHA(B, F_s, V_max, V_r, V_p):
    alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
    return alpha


def A(E, mu, ALPHA):
    a = 2 * PI * (2 * mu * E) ** 0.5 / ALPHA
    return a


def B(E, mu, V_p, V_r, ALPHA):
    b = 2 * PI * (2 * mu * (E - (V_p - V_r))) ** 0.5 / ALPHA
    return b


def D(bee, mu, ALPHA):
    d = 2 * PI * abs((2 * mu * bee - (ALPHA / 2) ** 2)) ** 0.5 / ALPHA
    return d


# Calculation of Transmission Probabilty of Eckart Potential
def T(a, b, d):
    T = (math.cosh(a + b) - math.cosh(a - b)) / (math.cosh(a + b) + math.cosh(d))
    return T


# Calculation of SINH function of Kappa
def S(V_max, E, TEMPERATURE):
    S = math.sinh(((V_max - E)) / (TEMPERATURE * k))
    return S


def sec_to_time(seconds):
    if seconds < 1:
        out = str(round(seconds,4)) + ' seconds'
        return out 
    else:
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        # suppose 365 days = 1 year, 30days = 1 month
        years = days // 365
        months = (days % 365) // 30
        days = (days % 365) % 30
        result = []
        count = 0
        if years:
            result.append(f"{years} years")
            count += 1
        if months and count < 2:
            result.append(f"{int(months)} months")
            count += 1
        if days and count < 2:
            result.append(f"{int(days)} days")
            count += 1
        if hours and count < 2:
            result.append(f"{int(hours)} hours")
            count += 1
        if minutes and count < 2:
            result.append(f"{int(minutes)} minutes")
            count += 1
        if seconds and count < 2:
            result.append(f"{int(seconds)} seconds")
            count += 1
        return " ".join(result)
    
def consumption_time(percent, k_r):
    # unit
    # k_rac [1/s]
    # consumption_time [s]
    consumption_ratio = percent / 100.0
    consumption_time = - math.log(1.0-consumption_ratio) / k_r
    return consumption_time

def pg_from_xyzfile(xyz_file):
    with open(xyz_file, "r") as file:
        file_content = file.read()
    first_line = file_content.strip().split('\n')[0]
    if len(first_line.split()) == 1:
        xyz_data = file_content.strip().split('\n')[2:]
    elif len(first_line.split()) == 4:
        xyz_data = file_content.strip().split('\n')
    syms = []
    poss = []
    for element in xyz_data:
        symbol = element.split()[0]
        syms.append(symbol)
        position = [float(flt) for flt in element.split()[1:]]
        poss.append(position)
    point_group = PointGroup(positions = poss, symbols = syms)
    return point_group.get_point_group()

input_file = sys.argv[1]
input = open(input_file, 'r')
input.seek(0)
input = list(input)

for line in input:
    if not line.strip() or line.split()[0] == '#':
        pass
    else:
        if 'Temperatures list' in line:
            Temp = [float(tmp) for tmp in line.split()[4:]]
            if len(Temp) == 1:
                pass
            elif Temp[-1] > Temp[-2]:
                pass                    
            elif Temp[-1] < Temp[-2] and len(Temp) == 3:
                Temp = [tmp for tmp in np.arange(Temp[0], Temp[1]+Temp[2], Temp[2])]
        elif 'Special Temperatures' in line:
            if len(line.split()) == 4:
                STemp = []
            elif line.split()[-1] == 'auto':
                STemp = [273.15, 293.15, 298.15]
            else:
                STemp = [float(tmp) for tmp in line.split()[4:]]
        elif 'RXN Symmetry Number' in line:
            Sym_num_in = line.split()[4]
            if Sym_num_in == 'auto':
                pass
            else:
                Sym_num = float(Sym_num_in)
        elif 'consumption percent' in line:
            consumption_percent = float(line.split()[3])

        elif 'Reactant Frequency file' in line:
            reactant_freq_path = line.split()[4]
        elif 'Product Frequency file' in line:
            product_freq_path = line.split()[4]
        elif 'TS Frequency file' in line:
            TS_freq_path = line.split()[4]
        elif 'Reactant Single Point file' in line:
            reactant_sp_path = line.split()[5]
        elif 'Product Single Point file' in line:
            product_sp_path = line.split()[5]
        elif 'TS Single Point file' in line:
            TS_sp_path = line.split()[5]
        elif 'Reactant XYZ file' in line:
            reactant_xyz_path = line.split()[4]
        elif 'Product XYZ file' in line:
            product_xyz_path = line.split()[4]
        elif 'TS XYZ file' in line:
            TS_xyz_path = line.split()[4]

Temp.extend(STemp)
Temp_list = sorted(Temp)

# Rxn symmetry
if pointgroup_install == 'YES' and Sym_num_in == 'auto':
    reactant_pg = pg_from_xyzfile(reactant_xyz_path) # Reactant point group
    reactant_rsn = rot_symmetry_number[reactant_pg] # Reactant rotational symmetry number
    TS_pg = pg_from_xyzfile(TS_xyz_path) # TS point group
    TS_rsn = rot_symmetry_number[TS_pg] # TS rotational symmetry number
    Sym_num = float(reactant_rsn) / float(TS_rsn)
elif pointgroup_install == 'NO' and Sym_num_in == 'auto':
    print("\n Fatal error : Invalid 'RXN Symmetry Number' in input file.\n")
    print(" RXN Symmetry Number 'auto' mode is depends on 'pointgroup' package")
    print(" Plz install 'pointgroup' package. ( pip install pointgroup )\n")
    sys.exit()

V_r = float(final_scf_energy(reactant_sp_path)) + float(zero_point_energy(reactant_freq_path))
V_p = float(final_scf_energy(product_sp_path)) + float(zero_point_energy(product_freq_path))
V_max = float(final_scf_energy(TS_sp_path)) + float(zero_point_energy(TS_freq_path))
G_r = float(final_scf_energy(reactant_sp_path)) + float(g_corr(reactant_freq_path))
G_p = float(final_scf_energy(product_sp_path)) + float(g_corr(product_freq_path))
G_ts = float(final_scf_energy(TS_sp_path)) + float(g_corr(TS_freq_path))

if V_r > V_p:
    E_o = V_r
else:
    E_o = V_p

# Scaling of Energies(define V_r == 0)
V_max = V_max - V_r
V_p = V_p - V_r
E_o = E_o - V_r
V_r = V_r - V_r
# V_p = 0.000/autokcal
# V_r = 0.00
# V_max = 9.80/autokcal
# E_o = 0.00/autokcal
y = (V_max - E_o) / 2.0
z = (V_max + E_o) / 2.0

# Specifing Parameters for the Eckart Potential
mu = float(reduced_mass(TS_freq_path)) * 1836
F_s = float(force_constant(TS_freq_path)) / 15.569141
bee = Bee(V_max, V_r, V_p)
alpha = ALPHA(bee, F_s, V_max, V_r, V_p)
d = D(bee, mu, alpha)


head ='TEMPERATURE      kappa_wigner     kappa_skodje     kappa_eckart     k_wiger             k_skodje              k_eckart             k_uncorrtd \n'
line = '------------------------------------------------------------------------------------------------------------------------------------------------'
calculation_info = f'\n# Level of thoery Info. # \nFreqency calculation  : {level_of_theory_freq(TS_freq_path)} \nSP calculation        : {level_of_theory_SP(TS_sp_path)}\n'

if Sym_num_in == 'auto':
    rxn_info = f'\n# RXN Info. # \nDetected Point group  : Reactant [ {reactant_pg} ] / TS [ {TS_pg} ] \nRXN Symmetry Number   : {Sym_num} \nRxn Barrier           : {(G_ts - G_r) * 2625.5e3}  [J/mol]\n\n'
    output = calculation_info + rxn_info + head + line + '\n'
else:
    thermodynamics = f'\n# Thermodynamics # \n\n RXN barrier  :  {(G_ts - G_r) * 2625.5e3}  [J/mol]\n\n' 
    output = calculation_info + thermodynamics + head + line + '\n'

for TEMPERATURE in Temp_list:
    raw = ''
    # Calculation of Uncorrected Gibbs Free Energy Barrier in kJ/mol and Rate
    delta_g_dagg = (G_ts - G_r) * 2625.5
    ln_uncorr_rate = math.log(BOLTZMANN_CONSTANT * TEMPERATURE / (PLANCK_CONSTANT)) - delta_g_dagg * 1000 / (
                GAS_CONSTANT * TEMPERATURE)
    uncorr_rate = "{:.4E}".format(math.exp(ln_uncorr_rate) * Sym_num)   ## uncorrected rate constant

    # Calculation of Wigner tunneling correction
    kappa_w = 1 + (((F_s / mu) ** 0.5 / (k * TEMPERATURE)) ** 2) / 24  ##   wigner correction
    
    # Calculation of Skodje tunneling correction
    alpha_s = (2 * PI / ((F_s / mu) ** 0.5))
    beta_s = (1 / (k * TEMPERATURE))
    if V_p < V_r:
        Vee = 0
    else:
        Vee = V_p - V_r
    if alpha_s > beta_s:
        kappa_s = (beta_s * PI / alpha_s) / (
                    math.sin(beta_s * PI / alpha_s)) + (beta_s / (alpha_s - beta_s) * math.exp(
            (beta_s - alpha_s) * (V_max - Vee)))
    else:
        kappa_s = (beta_s / (beta_s - alpha_s) * (math.exp((beta_s - alpha_s) * (V_max - Vee)) - 1))  ## Skodje correction

    # Calculation of Eckart tunneling correction using 10-point Gauss-Legendre Quadrature
    kappa = 1
    for i in range(10):
        a = A((x[i] * y + z), mu, alpha)
        b = B((x[i] * y + z), mu, V_p, V_r, alpha)
        kappa = (2 * y / (TEMPERATURE * k) * w[i] * S((V_max), (x[i] * y + z), TEMPERATURE) * T(a, b, d)) + kappa  ## Eckart correction

    # Calculation of Wigner Rate
    corr_rate_w = "{:.4E}".format(Sym_num * kappa_w * math.exp(ln_uncorr_rate))  ## wigner corrected rate constant

    # Calculation of Skodje Rate
    corr_rate_s = "{:.4E}".format(Sym_num * kappa_s * math.exp(ln_uncorr_rate))  ## Skodje corrected rate constant

    # Calculation of Eckart Rate
    corr_rate = "{:.4E}".format(Sym_num * kappa * math.exp(ln_uncorr_rate)) ## Eckart corrected rate constant

    raw = '           '.join(map(str, ['%.2f' %TEMPERATURE, '%.4f' %np.round(kappa_w,4), '%.4f' %np.round(kappa_s,4), '%.4f' %np.round(kappa,4), corr_rate_w,  corr_rate_s, corr_rate, uncorr_rate]))
    output = output + raw + '\n'

print(output)


with open(f"{input_file.split('.')[0]}.out", 'w') as file:
    file.write(output)