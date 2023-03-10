import numpy as np
import matplotlib.pyplot as plt
from Constants.MissionInputs import rho_0,rho_cruise,lbs_kg,V_cruise
from Constants.Aerodynamics import CL_Max_Clean
from Constants.Masses_Locations import m_mto,W_S_design

# Density
rho = 0 # just to define rho, this NEVER changes
density = 0 # 0 @ sea-level, 1 @ cruise, else @ loiter
if density == 0:
    rho = rho_0
elif density == 1:
    rho = rho_cruise
else:
    rho = 0.663838

if density == 0:
    rho = rho_0 # density at sea level
elif density == 1:
    rho = rho_cruise
else:
    rho = 0.663838

#----------------------------------------MANEUVER DIAGRAM DESIGN-------------------------------------------------------
#Max lift coefficient
if density == 0:
    Clmax = 1.45 #CL_max_landing
elif density == 1:
    Clmax = CL_Max_Clean
else:
    Clmax = CL_Max_Clean # add loiter Clmax to Constants file

print(Clmax)

#Load factor
m_mto_lbs = m_mto / lbs_kg #conversion between kg and lbs
n = 2.1 + (24000 / (m_mto_lbs + 10000))

if n < 2.5:
    n_max = 2.5
elif n > 3.8:
    n_max = 3.8
else:
    n_max = n
print(n_max)
n_min = -1

#Cruise speed EAS
V_C = V_cruise
V_C = V_C * np.sqrt(rho/rho_0)

#Stall speed EAS
V_S = np.sqrt((2 * 1 * W_S_design) / (rho * Clmax)) * np.sqrt(rho / rho_0)

#Dive speed EAS
if density == 0:
    h = 0
    dt = dt_takeoff
elif density == 1:
    h = h_cruise
    dt = dt_cruise
else:
    h = h_loiter
    dt = dt_loiter

a = ISA_calculator(h,dt)[3]

isa_cruise = ISA_calculator(h_cruise, dt_cruise)
# print("ISA CRUISE", isa_cruise)

M_C = V_C / a
M_D = M_C + 0.05

V_D1 = V_C/0.8
V_D2 = M_D * a
V_D = min(V_D1, V_D2)

# print("MACH NUMBER DIVE", V_D/a)