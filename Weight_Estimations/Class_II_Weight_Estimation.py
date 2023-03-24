from Constants.MissionInputs import ISA_calculator,h_cruise,lbs_kg,ft_m,g,FL_ft,V_cruise,dt_cruise,Aw,m_inches,rho_5000,V_approach
from Constants.Masses_Locations import W_P_design,m_f,m_oem,beta_s_land_fc
from Constants.FlightPerformance_Propulsion import eta_inverter,eta_EM,eta_wire,inverter_power_density, n_ult_pos
from Constants.Aerodynamics import CL_CD_DesCruise
from Constants.AircraftGeometry import bw, l_f,taperw, S_w, t_c_ratio_w, Sweep_quarterchordw, S_flap, d_f_outer
from Constants.Empennage_LandingGear import Sh,lh,A_h,bh,taperh, Sv, taperv, x_h, Av, tc_tail, Sweep_quarter_chord_H #, Sweep_quarter_chord_VT
from Constants.Aerodynamics import CL_MaxLand
import numpy as np



T_cruise, p_cruise, rho_cruise, a_cruise = ISA_calculator(h_cruise,0)
m_mto = 20000
no_fc = 11
m_tanks = 850
#all the equations from Raymer
#Calculated seperaterly:

#LH2 Storage/Fuel System
#Estimated at 1.4 grav index, add more detail later I would say
LH2_system_tank = m_tanks       #Must be done, thus change

#Fuel Cell
Fuel_Cell_Weight = no_fc * 80

#Wing
Wdg = m_mto/ lbs_kg
S_W = S_w * (1/ft_m)**2
Nz = n_ult_pos
t_cw = t_c_ratio_w
lambda_ = taperw
Lambda = Sweep_quarterchordw
S_csw = S_flap * (1/ft_m)**2
W_wing_lbs = 0.0051 * (Wdg*Nz)**0.557 * S_W**0.649 * Aw **0.5 * t_cw**(-0.4) * (1+lambda_)**0.1 * np.cos(Lambda)**(-1) * S_csw**0.1
W_wing = W_wing_lbs * lbs_kg

#Horizontal Tail
K_uht = 1.143
F_w = 1                     #fuselage width at horizontal tail intersection [ft] #TODO
B_h = bh * 1/ft_m           #horizontal tail span [ft]
L_t = 1                     #tail length; wing quarter-MAC to tail quarter-MAC [ft] #TODO
K_y = 1                     #aircraft pitching radius of gyration [ft] ( = 0.3Lt) #TODO
S_e = 1                     #elevator area [ft^2] #TODO
Sht = Sh * (1/ft_m)**2      #Horizontal Tail area [ft^2]
W_hor_tail_lbs = 0.0379 * K_uht *(1+ F_w/B_h)**(-0.25) * Wdg**0.639 * Nz**0.10 * Sh**0.75 * lh**(-1) * K_y* np.cos(Sweep_quarter_chord_H)**(-1) * A_h**0.166 * (1+S_e/Sht)**0.1
W_hor_tail = W_hor_tail_lbs * lbs_kg

# #Vertical Tail
H_t_H_v = 1.                # T tail or not, if T tail = 1 , if not = 0
S_vt = Sv/(ft_m)**2
K_z = 1                     #aircraft yawing radius of gyration, ft ( = L,)
W_ver_tail_lbs = 0.0026 * (1+ H_t_H_v)**0.225 * Wdg**0.556 * Nz **0.536 * L_t**(-0.5) * S_vt**0.5 * K_z**0.875 * np.cos(0)**(-1) * Av**0.35 * t_cw**(-0.5)
W_ver_tail = W_ver_tail_lbs *lbs_kg


#Fuselage
K_door = 1                  # 1.0 if no cargo door; = 1.06 if one side cargo door; = 1.12 if two side cargo doors; = 1.12 if aft clamshell door; = 1.25 if two side cargo doors and aft clamshell door
K_lg =  1                   # 1.12 if fuselage-mounted main landing gear;= 1.0 otherwise
L_f = l_f * (1/ft_m)
K_ws = 0.75 * ((1+2*taperw)/(1+taperw)) * (bw * np.tan(np.radians(Sweep_quarterchordw)/L_f))
S_f = 1                     # Fuselage wetted area [ft^2]
D_f = d_f_outer * (1/ft_m)   # Fuselage Depth/diameter [ft]
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
N_nw = 2                        # Number of nose wheels
W_noselg_lbs = 0.032 * K_np * W_l ** 0.646 * N_l ** 0.2 * L_n** 0.5 * N_nw**0.45
W_noselg = W_noselg_lbs * lbs_kg


# # Nacelle Group
power_req = m_mto * g / W_P_design
Eng_W_kg = 200*10**3 / 13 #Check Later
Engine_weight = power_req / eta_EM / Eng_W_kg + (power_req / eta_EM / eta_wire / eta_inverter) / inverter_power_density / 10**3
K_ng = 1.017                                #Pylon mounted nacelle or not, if not 1, if so 1.017
N_lt = 1                                    #Nacelle Length [ft]
N_w = 1                                     #Nacelle Width [ft]
W_ec = 2.331 * (Engine_weight/4 * (1/lbs_kg))**0.901 * 1.4         #Weight of engine and content [lbs]
N_en = 4                                    #Number of engines
S_n = 1                                     #Nacelle wetted area [ft^2]
W_nacelle_lbs = 0.6724 * K_ng * N_lt**0.10 * N_w**0.294 * Nz**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224
W_nacelle = W_nacelle_lbs * lbs_kg

# Engine Controls
L_ec = 1                                    #length from engine front to cockpit-total if multiengine [ft]
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

#
# #Hydrolics
# W_hydraulics_lbs = 0.008*W_dg                #CHANGE THIS LATER
# W_hydraulics = W_hydraulics_lbs * lbs_kg
#
# #Avionics
# #W_auv uninstalled avionics weight
# W_uav = 800
#
# W_av_lbs = 2.117*W_uav**(0.933)
# W_av = W_av_lbs * lbs_kg
#
# #Electrical
# W_elec_lbs = 12.57*(W_av_lbs)**(0.51)
# W_elec = W_elec_lbs * lbs_kg
#
# #Airconditioning and Anti-Icing
# #N_p personnel
# N_p = 3
# #Mach Number
# Mach = V_cruise/a_cruise
#
# W_aircon_ice_lbs = 0.265*W_dg**(0.52)*N_p**(0.68)*W_av_lbs**(0.17)*Mach**(0.08)
# W_aircon_ice = W_aircon_ice_lbs * lbs_kg
#
# #Furnishing
# W_furnish_lbs = 0.0582*W_dg-65
# W_furnish = W_furnish_lbs * lbs_kg
#
# #Paint
# W_paint_lbs = 0.0045 * W_dg * 0
# W_paint = W_paint_lbs * lbs_kg

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
# print('Weight Flight Controls =', W_flight_controls)
# print('Weight Hydraulics =', W_hydraulics)
# print('Weight Avionics =', W_av)
# print('Weight Electronics =', W_elec)
# print('Weight Aircon and De-Icing', W_aircon_ice)
# print('Weight Furnishing =', W_furnish)
# print("Weight Installed Engines =", W_installed_eng)
# print("Weight Paint =", W_paint)
# all_masses = [LH2_system_tank, W_nacelle, Fuel_Cell_Weight, W_wing, W_hor_tail, W_ver_tail, W_fus, W_land_main, W_land_nose, W_flight_controls, W_hydraulics, W_av, W_elec, W_aircon_ice, W_furnish, W_installed_eng, W_paint]
# print("Fuel Mass", m_f)
# print("MTOM", m_mto)
# print('Class II Weight Estimation =', LH2_system_tank+Engine_weight+Fuel_Cell_Weight+W_wing+W_hor_tail+W_ver_tail+W_fus+W_land_main+W_land_nose+W_flight_controls+W_hydraulics+W_av+W_elec+W_aircon_ice+W_furnish+W_installed_eng+W_paint)
# print('Class I OEM =', m_oem)
# print(sum(all_masses))
# print(all_masses)