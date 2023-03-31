from Constants.FlightPerformance_Propulsion import eta_EM, eta_wire, eta_fuelcell, eta_inverter, Power_tot, Voltage_cell
from Constants.MissionInputs import PAX, V_cruise, rho_cruise, V_approach, rho_5000, rho_0, g
from Constants.Masses_Locations import m_mto
from Constants.AircraftGeometry import S_w
from Constants.Aerodynamics import CL_MaxLand, CL_MaxTakeOff
import numpy as np
# Preliminary Calculations
eta_tot = eta_fuelcell * eta_inverter * eta_wire * eta_EM
V_TakeOff_0 = np.sqrt(m_mto*g/(0.5*rho_0*S_w*CL_MaxTakeOff))
V_TakeOff_5000 = np.sqrt(m_mto*g/(0.5*rho_5000*S_w*CL_MaxTakeOff))
V_Land_0 = 1.15*np.sqrt(m_mto*g/(0.5*rho_0*S_w*CL_MaxLand))
V_Land_5000 = 1.15*np.sqrt(m_mto*g/(0.5*rho_5000*S_w*CL_MaxLand))

# Constants
MM_O2 = 32
n_electrons = 4
Faraday = 96485
g_O2 = 0.23133
lambda_O2 = 2.8

#Cabin Pressurization  #FAA: 0.25kg/min/pax of air  https://www.faa.gov/newsroom/cabin-air-quality#:~:text=FAA%20regulations%20require%20airliners'%20ventilation,consistent%20with%20other%20public%20environments.

mdot_cabin = 0.25 * (PAX+3)  / 60

# Mass flow of Air
mdot_O2 = MM_O2 * Power_tot/ n_electrons/ Faraday/ Voltage_cell # todo: Take off/cruise differences
mdot_air = mdot_O2 * lambda_O2 /g_O2/ eta_tot/1000

# Mass flow for cooling
mdot_cooling = 0

# Total Mass flow
mdot_total = mdot_air + mdot_cabin + mdot_cooling

# Velocity - free stream velocity as it in front of fuselage
prodruding_ratio = 1
V_boundary = [V_TakeOff_0, V_TakeOff_5000, V_cruise, V_approach, V_Land_0, V_Land_5000]  # todo: check for take off etc
rho_phases = [rho_0, rho_5000, rho_cruise, rho_0, rho_5000, rho_5000]
FlightPhases = ["Take off, SL","Take off, 5000", "Cruise",  "Approach, SL", "Landing, SL", "Landing, 5000"]
# IFF there is a protrusion, Use: (7/8)*prodruding_ratio**(8/7)*V_cruise

# Area of inlet
A_inlet = []
for i in range(len(V_boundary)):
    A = mdot_total / (V_boundary[i] * rho_phases[i])
    A_inlet.append(A)

A_inlet_critical = max(A_inlet)
FlightCondition = FlightPhases[A_inlet.index(max(A_inlet))]

# PRINT STATEMENTS
print("Total mass flow", mdot_total, "kg/s")

print("Inlet area options --------- ")
for i in range(len(V_boundary)):
    print(A_inlet[i], "m^2 during", FlightPhases[i])
print("----------------")

print("Inlet area maximum required", A_inlet_critical, "m^2", "during", FlightCondition)
