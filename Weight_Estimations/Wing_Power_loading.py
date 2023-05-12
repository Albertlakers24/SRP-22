import numpy as np
from matplotlib import pyplot as plt
from Constants.MissionInputs import ISA_calculator,s_takeoff,s_landing,h_cruise,dt_cruise,rho_0,takeoff_critical,dt_takeoff,Aw,V_approach,ft_m,V_cruise,g
from Constants.FlightPerformance_Propulsion import eta_prop
from Constants.Aerodynamics import CL_MaxLand,CL_MaxTakeOff,CD_DesTakeOff,CD0_40,CD0_CR
#Constants
#Constants
rho_1524= ISA_calculator(takeoff_critical,dt_takeoff)[2]              #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise = ISA_calculator(h_cruise,dt_cruise)[2]
W_S = np.arange(1,6000,1)
s_takeoff_1524 = s_takeoff
s_landing_1524 = s_landing
W_S = np.arange(1,6000,0.01)
Psi = 0.0075                        #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                          #span efficiency factor (value based on Roelof reader p.46)
Cfe = 0.0030                        #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                        #(6.0-6.2) wetted area ratios -> depending on airframe structure
##Cdo calculations
e = 1/(np.pi*Aw*Psi+(1/phi))
#ROC and beta estimates
ROC = 6                        #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_OEI = 0.5
ROC_V = 0.032#0.0032
ROC_V_OEI = 0.03                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach_stall = V_approach /1.23  #CS 25 requirement of V_stall_land = V_approach / 1.23
beta_V_app_fc = 0.985
beta_s_land_fc = 0.95
beta_cruise_fc  = 0.99
beta_ROC_fc = 0.99
beta_cv = 1
beta_s_to = 1
beta_em = 1
C_LFL = 0.45                        #Landing field length coefficient s^2/m
alpha_p_em = 1
CL_2 = (1/1.13)**2 * CL_MaxTakeOff
k_t = 0.85
h2 = 50 * ft_m
V_s0 = 31
def oswald_efficiency(flap_deflection):
    delta_e = 0.0026 * flap_deflection
    e_new = e + delta_e
    return e_new
def Power_lapse(rho,rho_0):
    alpha_p_ice = (rho/rho_0) ** (3/4)
    return alpha_p_ice
def V_approach_constraint(density,V_approach_stall,CL_MaxLand,beta):
    W_S = 1/2 * density * V_approach_stall **2 * CL_MaxLand * 1/beta
    return W_S
def s_land_constraint(s_land,C_LFL,density,CL_MaxLand,beta):
    W_S = (s_land / C_LFL) * (density * CL_MaxLand / 2) * 1/beta
    return W_S
def cruise_contraint(eta_prop,alpha_p,Cd0,density,V,W_S,Aw,e,beta):
    W_P = eta_prop * (alpha_p/beta) *((Cd0*1/2*density*V**3)/(beta*W_S)+ (beta*W_S)/(np.pi*Aw*e*1/2*density*V))**(-1)
    return W_P
def roc_constraint(eta,alpha_p,ROC,Cd0,density,Aw,e,W_S,beta,N_e,y):
    if y == 1: #One engine inoperative case
        W_P = ((N_e - 1)/N_e)*eta * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*Aw*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    else:
        W_P = eta * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*Aw*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    return W_P
def climb_gradient_constraint(eta,alpha_p,ROC_V,CD,CL,density,W_S,beta,N_e,y):
    if y == 1:
        W_P = ((N_e - 1)/N_e)* eta * (beta / alpha_p) * (1 / (ROC_V + (CD / CL))) * np.sqrt((density / 2) * ((CL) / (beta * W_S)))
    else:
        W_P = eta_prop * (beta/alpha_p) * (1/(ROC_V + (CD/CL))) * np.sqrt((density/2)*((CL)/(beta*W_S)))
    return W_P
def takeoff_constraint(alpha_p,L_to,density,h_2, k_t,N_e,y):
    if y == 1:
        W_P = alpha_p  * (1.15 * np.sqrt((N_e/(N_e-1))*(W_S/(L_to * k_t * density * g * np.pi * Aw * e))) + (N_e/(N_e-1))*(4*h_2/L_to)) **(-1) * np.sqrt((CL_2/W_S) *(density)/2)
    else:
        W_P = alpha_p * ((1.15 * np.sqrt((W_S / (L_to * k_t * density * g * np.pi * Aw * e)))) + (4 * h_2 / L_to)) ** (-1) * np.sqrt((CL_2 / W_S) * ((density) / 2))
    return W_P
e_to = oswald_efficiency(40)
W_S_approach = V_approach_constraint(rho_1524, V_approach_stall, CL_MaxLand, beta_V_app_fc)
W_S_land = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_MaxLand, beta_s_land_fc)
W_P_cruise = cruise_contraint(eta_prop,alpha_p_em, CD0_CR, rho_cruise, V_cruise, W_S, Aw, e,beta_cruise_fc)
W_P_ROC = roc_constraint(eta_prop, alpha_p_em, ROC, CD0_40, rho_1524, Aw, e_to, W_S, beta_ROC_fc, 4, 2)
W_P_ROC_OEI = roc_constraint(eta_prop, alpha_p_em, ROC_OEI, CD0_CR, rho_cruise, Aw, e, W_S, beta_cruise_fc, 4, 1)
W_P_CV = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V, CD_DesTakeOff, CL_MaxTakeOff, rho_1524, W_S,beta_em, 4, 2)
W_P_CV_OEI = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V_OEI, CD_DesTakeOff, CL_MaxTakeOff, rho_1524,W_S, beta_cv, 4, 1)
W_P_TOP = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 4, 2)
W_P_TOP_OEI = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 4, 1)
#DESIGN POINTS OTHER AIRCRAFT
#ATR 42
power_42 = 2 * 1610710
weight_42 = 18600 * g
surface_42 = 54.5
WP_42 = weight_42 / power_42
WS_42 = weight_42 / surface_42
#ATR 72
power_72 = 1610710 * 2
weight_72 = 22800 * g
surface_72 = 61
WP_72 = weight_72 / power_72
WS_72 = weight_72 / surface_72
#DASH 8
power_dash = 1864250 * 2
weight_dash = 19505 * g
surface_dash = 56.3
WP_dash = weight_dash / power_dash
WS_dash = weight_dash / surface_dash
print(WS_dash, WP_dash)

plt.vlines(W_S_approach,0,100,'b',label="Approach Speed Constraint(141kts)")
plt.plot(W_S,W_P_TOP,'g',label = "Takeoff Constraint, CL = 1.9")
plt.plot(W_S,W_P_TOP_OEI,'r',label = "Takeoff Constraint (OEI), CL = 1.9")
plt.vlines(W_S_land,0,100,'c',label ="Landing Constraint, CL = 2.2")
plt.plot(W_S,W_P_cruise,'m',label = "Cruise Constraint (275kts)")
plt.plot(W_S,W_P_ROC,'y',label = "Rate of Climb Constraint, ROC = 6")
plt.plot(W_S,W_P_ROC_OEI,'orange',label = "Rate of Climb Constraint (OEI), ROC = 0.5")
plt.plot(W_S,W_P_CV,'k',label = "Climb Gradient Constraint, $\gamma=0.032$")
plt.plot(W_S,W_P_CV_OEI,'indigo',label = "Climb Gradient Constraint (OEI), $\gamma=0.03$")
plt.plot(WS_42, WP_42, "*", label="Design Point ATR 42", markersize = 7)
plt.plot(WS_72, WP_72, "^", label="Design Point ATR 72", markersize = 7)
plt.plot(WS_dash, WP_dash, "s", label="Design Point Dash 8 q300", markersize = 7)
plt.xlim(0,6000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,6001,500),fontsize=12)
plt.yticks(np.arange(0,0.5,0.05),fontsize=12)
plt.xlabel("W/S (N/$m^2$)",fontsize=16)
plt.ylabel("W/P (N/W)",fontsize=16)
plt.fill_between(W_S, W_P_TOP_OEI, 1, color="red", alpha=0.1)
plt.axvspan(W_S_land, 6000, color="red", alpha=0.1)
plt.fill_between(W_S, W_P_cruise, 1, color="red", alpha=0.1)
plt.fill_between(W_S, W_P_CV_OEI, 1, color="red", alpha=0.1)
plt.fill_between(W_S, W_P_ROC, 1, color="red", alpha=0.1)
plt.plot(W_S[281960], W_P_TOP_OEI[281960], 'o',label = "Design Point SRP-22", markersize = 7)
plt.legend(loc = "upper right")
plt.grid()
plt.show()