from Constants.MissionInputs import ISA_calculator,h_cruise,lbs_kg,ft_m,g,Aw,m_inches,rho_5000,V_approach
from Constants.Masses_Locations import W_P_design,m_f,m_pldes, m_mto, m_oem
from Constants.FlightPerformance_Propulsion import eta_inverter,eta_EM,eta_wire,inverter_power_density, n_ult_pos
from Constants.AircraftGeometry import bw, l_f,taperw, S_w, t_c_ratio_w, Sweep_quarterchordw, S_flap, d_f_outer
from Constants.Empennage_LandingGear import Sh,A_h,bh, Sv, lh, Av, Sweep_quarter_chord_HT, Sweep_halfchord_VT, lv
from Constants.Aerodynamics import CL_MaxLand
import numpy as np



T_cruise, p_cruise, rho_cruise, a_cruise = ISA_calculator(h_cruise,0)
no_fc = 11
m_tanks = m_f/0.2
#all the equations from Raymer
#Calculated seperaterly:

#LH2 Storage/Fuel System
#Estimated at 1.4 grav index, add more detail later I would say
LH2_system_tank = m_tanks      #Must be done, thus change

#Fuel Cell
Fuel_Cell_Weight = no_fc * 80

#Wing
Wdg = m_mto * (1 / lbs_kg)
S_W = S_w * (1/ft_m)**2
Nz = n_ult_pos
t_cw = t_c_ratio_w
lambda_ = taperw
Lambda = Sweep_quarterchordw
S_csw = (S_flap + 5) * (1/ft_m)**2
W_wing_lbs = 0.0051 * (Wdg*Nz)**0.557 * S_W**0.649 * Aw**0.5 * t_cw**(-0.4) * (1+lambda_)**0.1 * np.cos(np.radians(Lambda))**(-1) * S_csw**0.1
W_wing = W_wing_lbs * lbs_kg

#Horizontal Tail
K_uht = 1.143
F_w = 0 * (1/ft_m)                     #fuselage width at horizontal tail intersection [ft]
B_h = bh * (1/ft_m)                    #horizontal tail span [ft]
L_t = lh * (1/ft_m)                    #tail length; wing quarter-MAC to tail quarter-MAC [ft]
K_y = 0.3 * L_t                        #aircraft pitching radius of gyration [ft] ( = 0.3Lt)
S_e = 1.35 * (1/ft_m)**2               #elevator area [ft^2]
Sht = Sh * (1/ft_m)**2                 #Horizontal Tail area [ft^2]
W_hor_tail_lbs = 0.0379 * K_uht * (1+F_w/B_h)**(-0.25) * Wdg**(0.639) * Nz**(0.10) * Sht**(0.75) * L_t**(-1.0) * K_y**(0.704) * np.cos(Sweep_quarter_chord_HT)**(-1.0) * A_h**(0.166) * (1+S_e/Sht)**(0.1)
W_hor_tail = W_hor_tail_lbs * lbs_kg

# #Vertical Tail
H_t_H_v = 1                # T tail or not, if T tail = 1 , if not = 0
S_vt = Sv * (1/ft_m)**2
L_v = lv * (1/ft_m)
K_z = L_v                     #aircraft yawing radius of gyration, ft ( = Lt)
t_c_v = 0.18
W_ver_tail_lbs = 0.0026 * (1+ H_t_H_v)**0.225 * Wdg**0.556 * Nz **0.536 * L_v**(-0.5) * S_vt**0.5 * K_z**0.875 * np.cos(Sweep_halfchord_VT)**(-1) * Av**0.35 * t_c_v**(-0.5)
W_ver_tail = W_ver_tail_lbs * lbs_kg

#Fuselage
K_door = 1.06               # 1.0 if no cargo door; = 1.06 if one side cargo door; = 1.12 if two side cargo doors; = 1.12 if aft clamshell door; = 1.25 if two side cargo doors and aft clamshell door
K_lg = 1                    # 1.12 if fuselage-mounted main landing gear;= 1.0 otherwise
L_f = l_f * (1/ft_m)
K_ws = 0.75 * ((1+2*taperw)/(1+taperw)) * (bw * np.tan(np.radians(Sweep_quarterchordw)/L_f))
S_f = 105.5 * (1/ft_m)**2                 # Fuselage wetted area [ft^2]
D_f = d_f_outer * (1/ft_m)                # Fuselage Depth/diameter [ft]
L_D = L_f / D_f
W_fus_lbs = 0.3280 * K_door * K_lg * (Wdg*Nz)**0.5 * L_f ** 0.25 * S_f **0.302 * (1 + K_ws)** 0.04 * (L_D)**0.10
W_fus = W_fus_lbs * lbs_kg

# # #Main Landing Gear
K_mp = 1.0                      #Kneeling gear or not, 1 if not , if so 1.126
W_l = Wdg                       #landing design gross weight, lb
N_gear = (CL_MaxLand*1/2 *rho_5000 *V_approach**2 * S_w)/(m_mto* g) #Landing load factor
N_l = N_gear * 1.5              #ultimate landing load factor; = Ngear X 1.5
L_m = 3.5 * m_inches            #length of main landing gear [inches] 1.8m in belly/ 3.5m in nacelle
N_mw = 4                        #length of main landing gear
N_mss = 2                       #number of main wheels
V_stall = V_approach/1.23       #V_stall at landing
W_mainlg_lbs = 0.0106 * K_mp * W_l ** 0.888 * N_l ** 0.25 * L_m ** 0.4 * N_mw **0.321 * N_mss ** (-0.5) * V_stall ** 0.1
W_mainlg = W_mainlg_lbs * lbs_kg

# #Nose Landing Gear
K_np = 1.15                     #Kneeling gear or not, 1 if not , if so 1.15
L_n = 1.2 * m_inches            #Nose gear length [inches]  [1.2m]
N_nw = 2                        #Number of nose wheels
W_noselg_lbs = 0.032 * K_np * W_l ** 0.646 * N_l ** 0.2 * L_n** 0.5 * N_nw**0.45
W_noselg = W_noselg_lbs * lbs_kg

# # Nacelle Group
power_req = m_mto * g / W_P_design
Eng_W_kg = 5.5 * 10**3 #Check Later
Engine_weight = power_req / eta_EM / Eng_W_kg + (power_req / eta_EM / eta_wire / eta_inverter) / inverter_power_density / 10**3
K_ng = 1.017                                #Pylon mounted nacelle or not, if not 1, if so 1.017
N_lt = 5.14 * (1/ft_m)                       #Nacelle Length [ft]
N_w = 1 * (1/ft_m)                           #Nacelle Width [ft]
K_p = 1.4
W_ec = 2.331 * (Engine_weight/4 * (1/lbs_kg))**0.901 * K_p         #Weight of engine and content [lbs]
N_en = 4                                    #Number of engines
S_n = 22.76 * (1/ft_m)**2                                     #Nacelle wetted area [ft^2]
W_nacelle_lbs = 0.6724 * K_ng * N_lt**0.10 * N_w**0.294 * Nz**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224
W_nacelle = W_nacelle_lbs * lbs_kg

# Engine Controls
L_ec = 14.3 * (1/ft_m) * 4 * 1.5                         #length from engine front to cockpit-total if multiengine [ft]
W_engine_control_lbs = 5 * N_en + 0.8 * L_ec
W_engine_control = W_engine_control_lbs * lbs_kg

# Starter/ Pneumatics
W_en = Engine_weight/4 * (1/lbs_kg)
W_starter_lbs = 49.19 * ((N_en * W_en)/1000)**0.541
W_starter = W_starter_lbs * lbs_kg

#Flight Controls
N_f = 4
N_m = 2
S_cs = 1
I_y =  1                #yawing moment of inertia
W_flight_controls_lbs = 145.9 * N_f**0.554 * (1+N_m/N_f)**(-1) * S_cs **0.20 * (I_y * 10**(-6))**0.07
W_flight_controls = W_flight_controls_lbs * lbs_kg

#Instrument
Kr = 1                  # Reciprocating engine or not
Nc = 3                  # Number of crew
W_instrument_lbs = 4.509 * Kr * K_p * Nc**0.541 * N_en * (L_f + bw)**0.5
W_instrument = W_instrument_lbs * lbs_kg

#Hydraulics
W_hydraulics_lbs = 0.2673 * N_f * (l_f + bw)**0.937
W_hydraulics = W_hydraulics_lbs * lbs_kg

#Electrical
R_kva = 60                         #system electrical rating, kv · A (typically 40-60 for transports)
L_a = 100 * (1/ft_m)               #electrical routing distance, generators to avionics to cockpit [ft]
N_gen = N_en                       #Number of Generator (typically the same as N_en)
W_electrical_lbs = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10
W_electrical = W_electrical_lbs * lbs_kg

#Avionics
W_uav = 1400                    #Typically around 800- 1400 lbs
W_avionics_lbs = 1.73 * W_uav**0.983
W_avionics = W_avionics_lbs * lbs_kg

#Furnishing
W_c = m_pldes * (1/lbs_kg)
W_furnishing_lbs = 0.0577*Nc**0.1 * W_c**0.393 * S_f**0.75
W_furnishing = W_furnishing_lbs * lbs_kg

#Air conditioning
N_p = 51
V_pr = 88.5 * (1/ft_m)**3                    #Volume of presurried section [ft^3]
W_airconditioning_lbs = 62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10
W_airconditioning = W_airconditioning_lbs * lbs_kg

#Anti-Icing
W_antiice_lbs = 0.002 * Wdg
W_antiice = W_antiice_lbs * lbs_kg

#Handling Gear
W_handling_gear_lbs = 3.0 * 10**(-4) * Wdg
W_handling_gear = W_handling_gear_lbs * lbs_kg
#Results:

#Things to be reconsidered later before iteration etc:

#Fuselage pressurized volume
#Landing gear load
#Both landing gear length

print('Results Start Here:')
print('Weight Storage System LH2 =', LH2_system_tank)
print('Weight Fuel Cell = ', Fuel_Cell_Weight)
print('Weight Wing =', W_wing)
print('Weight Horizontal Tail = ', W_hor_tail)
print('Weight Vertical Tail = ', W_ver_tail)
print('Weight Fuselage =', W_fus)
print('Weight Main Landing Gear =', W_mainlg)
print('Weight Nose Landing Gear =', W_noselg)
print('Weight Engine Group =', W_nacelle)
print('Weight Engine Controls =', W_engine_control)
print('Weight Flight Control =', W_flight_controls)
print('Weight Instrument =', W_instrument)
print('Weight Hydraulics =', W_hydraulics)
print('Weight Electrical', W_electrical)
print('Weight Avionics', W_avionics)
print('Weight Furnishing =', W_furnishing)
print("Weight Air Conditioning =", W_airconditioning)
print("Weight Anti-Icing =", W_antiice)
print("Weight Handling Gear =", W_handling_gear)
all_masses = [LH2_system_tank, Fuel_Cell_Weight,W_wing,W_hor_tail,W_ver_tail,W_fus,W_mainlg,W_noselg,W_nacelle,W_engine_control,W_starter,W_flight_controls,W_instrument, W_hydraulics,W_electrical,W_avionics,W_furnishing,W_airconditioning,W_antiice,W_handling_gear]
print("Fuel Mass", m_f)
print("MTOM", m_mto)
print('Class II Weight Estimation =', sum(all_masses))
print('Class I OEM =', m_oem)