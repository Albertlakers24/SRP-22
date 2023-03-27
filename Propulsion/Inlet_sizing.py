from Constants.FlightPerformance_Propulsion import eta_EM, eta_wire, eta_fuelcell, eta_inverter, Power_tot, Voltage_cell
from Constants.MissionInputs import PAX, V_cruise, rho_cruise

# Imports
eta_tot = eta_fuelcell * eta_inverter * eta_wire * eta_EM

# Constants
MM_O2 = 32
n_electrons = 4
Faraday = 96485
g_O2 = 0.23133
lambda_O2 = 2.8

#Cabin Pressurization  #FAA: 0.25kg/min/pax of air  https://www.faa.gov/newsroom/cabin-air-quality#:~:text=FAA%20regulations%20require%20airliners'%20ventilation,consistent%20with%20other%20public%20environments.

mdot_cabin = 0.25 * (PAX+3)  / 60

# Mass flow of Air
mdot_O2 = MM_O2 * Power_tot/ n_electrons/ Faraday/ Voltage_cell
mdot_air = mdot_O2 * lambda_O2 /g_O2/ eta_tot/1000

# Mass flow for cooling
mdot_cooling = 0

# Total Mass flow
mdot_total = mdot_air + mdot_cabin + mdot_cooling

# Area of inlet - free stream velocity as it in front of fuselage
prodruding_ratio = 1
V_boundary =  V_cruise  # IFF there is a protrusion, Use: (7/8)*prodruding_ratio**(8/7)*V_cruise

A_inlet = mdot_total / (V_cruise * rho_cruise)

# PRINT STATEMENTS
print("Total mass flow", mdot_total, "kg/s")
print("Inlet area", A_inlet, "m^2")

