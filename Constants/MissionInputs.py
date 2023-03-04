import numpy as np
#Units Conversion
kts_m_s = 0.514444                  #knots to m/s
nmi_mile = 1.151                    #nmi to miles
nmi_m = 1852                        #nmi to m
min_s = 60                          #minutes to seconds
ft_m = 0.3048                       #ft to m
km_m = 1000                         #km to m
FL_ft = 100                         #FL to ft
lbs_kg = 0.453592                   #lbs to kg
kW_hp = 1.341                       #kW to hp
#General Constants
g = 9.80665                         #gravity
T_0 = 288.15                        #ISA Temperature
p_0 = 101325                        #ISA Pressure
rho_0 = 1.225                       #ISA Density
lapse_rate = -6.5/km_m              #Troposphere lapse rate
Molar_mass_air = 0.0289644          #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052     #J·kg⁻¹·K⁻¹
universal_gas_constant = 8.31432    #N m kmol⁻¹ K⁻¹
gamma = 1.4                         #Specific heat ratio of gas
# e_lh2 = 120*10**6                   #J/kg Specific Energy Liquid Hydrogen
#Requirment Constants
V_cruise = 275 * kts_m_s            #Cruise Velocity Requirement
h_cruise = 280 * FL_ft * ft_m       #Cruise height Requirement
h_loiter = 195 * FL_ft * ft_m       #Loiter height (requirement)
V_approach = 141 * kts_m_s          #Approach speed Requirement
takeoff_critical = 5000 * ft_m      #Takeoff at 5000ft above mean sea level
landing_critical = 5000 * ft_m      #Landing at 5000ft above mean sea level
dt_land = 10                        #Offset temperature to ISA at landing
dt_takeoff = 10                     #Offset temperature to ISA at takeoff
dt_cruise = 0                       #Offset temperature to ISA at cruise
dt_loiter = 0                       #Offset temperature to ISA at loiter
PAX = 48                            #Passenger Requirement
E = 30 * 60                         #Loiter endurance in seconds
s_takeoff = 4500 * ft_m             #Takeoff Distance
s_landing = 4500 * ft_m             #Landing Distance
m_pax = 200*lbs_kg*PAX              #Fat American passenger mass
m_pax_baggage = 40*lbs_kg*PAX       #Baggage mass for passenger
vol_pax = 0.14                      #Volume per passenger
vol_pax_baggage = 0.14 * PAX        #Total volume of baggage
m_crew = 3 * 190 * lbs_kg           #Crew mass
m_crew_baggage = 3 * 30 * lbs_kg    #Baggage mass for crew
R_norm = 1000 * nmi_m               #Design Range
R_div = 100 * nmi_m                 #Divergence Range
t_loiter = 30 * min_s               #Loiter Endurance
f_con = 5/100                       #Contingency fuel percentage
A = 12                            #Aspect Ratio (ONLY VALUE THAT COULD BE ITERATED)

''' THESE ARE IN PROPULSION SYSTEMS '''
#Efficiency Constants
# eta_prop = 0.85                     #Propeller efficiency
# eta_EM = 0.95                       #Electric motor efficiency
# eta_wire = 0.97                     #Wire efficiency
# eta_inverter = 0.995                #Inverter efficiency
# eta_fuelcell = 0.60                 #Fuel cell efficiency
# fc_power_density = 3                #kW/kg
# inverter_power_density = 30         #kW/kg
# em_power_density = 8                #kW/kg
# H2_power_density = 120              #MJ/kg

#Aerodynamic Constants
''' PLS ADD THE NEXT FEW COMMENTED TO WEIGHT ESTIMATION FILES SEPARATELY --- THESE ARE NOT MISSION INPUTS 
                    AND ARE ONLY USED THERE AS FAR AS I (MEGHA) KNOW  '''
# Psi = 0.0075                        #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
# phi = 0.97                          #span efficiency factor (value based on Roelof reader p.46)
# Cfe = 0.0030                        #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
# Swet_S = 6.1                        #(6.0-6.2) wetted area ratios -> depending on airframe structure

''' THESE ARE IN AERODYNAMICS '''
# CL_max_takeoff = 2.4                 # -
# CL_max_cruise = 1.5                 # -
# CL_max_landing = 2.4                # -
# CL_max_loiter = 1.5                 # -
''' THESE ARE IN MASSES_LOCATIONS '''
#Calculation Constants
# OEW_cg = 11.3                       #m

#General Functions
def ISA_calculator(h,dt):
    T = T_0 + lapse_rate * h + dt
    p = p_0 * (((T-dt) / T_0) ** ((-g) / (specific_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant*T)
    a = np.sqrt(gamma * T * specific_gas_constant)
    return T, p, rho, a

T_5000, p_5000, rho_5000, a_5000 = ISA_calculator(takeoff_critical, 10)
T_cruise, p_cruise, rho_cruise, a_cruise = ISA_calculator(h_cruise, 0)

kv_cruise = 3.046*10**(-5)          # Kinematic viscosity -- from ISA calculator website - correct it if its wrong
kv_TO_Land = 1.461 * 10**(-5)       # Kinematic viscosity

def Calculate_Reynolds_number(V, l, kv):
    Re = V*l/kv # kv = kinematic viscosity
    return Re