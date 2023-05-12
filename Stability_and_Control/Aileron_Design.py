import numpy as np

from Constants.MissionInputs import V_approach, g, ISA_calculator, takeoff_critical, dt_takeoff, rho_5000, V_cruise,a_cruise
from Constants.AircraftGeometry import S_w,bw,c_rw,c_tw,taperw
from Aerodynamics.AeroFunctions import Calculate_wingsweep
from Constants.Masses_Locations import m_mto
import pandas as pd
aileron_efficiency = 6.7/2.6 * 0.2
Ca_c = 0.3
aileron_percent = 0.535
aileron_lm_choice = aileron_percent * bw/2
V_approach_stall = V_approach /1.23  #CS 25 requirement of V_stall_land = V_approach / 1.23
V_stall = np.sqrt(m_mto * g / (ISA_calculator(takeoff_critical,dt_takeoff)[2]*0.5*2.4*S_w))
CL_alpha_AF = 6.86861805
Cd0_AF_ALpha0 = 0.0046
column_name = ["ld","l1","Ca_c","dA","dA_upper","dA_lower","dt"]
aileron_DP = pd.DataFrame(columns = column_name)
def yield_func(low,high,step):
    i = low
    while i <= high:
        yield i
        i+=step
def trapezoid_area(h,top,base):
    area = (top+base)/2 * h
    return area

def l_top(alpha,beta,height,base):
    x1 = height/ np.tan(alpha)
    x2 = height/ np.tan(beta)
    top = base - x1 - x2
    return top

def roll_control_derivative(innerboard,outerboard):
    outerboard_integral = c_rw/2*outerboard**2 - 2/3*(c_rw-c_tw)/(bw)*outerboard**3
    innerboard_integral = c_rw/2*innerboard**2 - 2*(c_rw-c_tw)/(3*bw)*innerboard**3
    Cl_dA = (2 * CL_alpha_AF * aileron_efficiency)/(S_w*bw) * (outerboard_integral - innerboard_integral)
    return Cl_dA

def roll_damping_coefficient(innerboard,outerboard):
    outerboard_integral = c_rw/3*outerboard**3 - (c_rw-c_tw)/(2*bw)*outerboard**4
    innerboard_integral =c_rw/3*innerboard**3 - (c_rw-c_tw)/(2*bw)*innerboard**4
    Cl_p = (-4 * (CL_alpha_AF + Cd0_AF_ALpha0))/(S_w*bw**2) * (outerboard_integral - innerboard_integral)
    return Cl_p
def steady_state_roll(Cl_dA,Cl_p,dA,V):
    P = -Cl_dA/Cl_p * np.radians(dA) * (2*V/bw)
    return P

def time_req(angle,P):
    dt = np.radians(angle)/P
    return dt

LE_hinge_line_angle_deg = Calculate_wingsweep(0,0.25,0.15,taperw)   #front spar position at 15% chord
TE_hinge_line_angle_deg = Calculate_wingsweep(0,0.25,0.60,taperw)    #rear spar position at 60% chord
LE_angle_deg = Calculate_wingsweep(0,0.25,0,taperw)                 #Leading edge angle
TE_angle_deg = Calculate_wingsweep(0,0.25,1,taperw)                 #Trailing edge angle
alpha_angle = np.radians(90 + TE_angle_deg)
beta_angle = np.radians(90 - LE_angle_deg)

l1_sequence = yield_func(aileron_percent*bw/2,0.9*bw/2,0.05)
for l1 in l1_sequence:
    dA_sequence = yield_func(0,23,0.1)
    for dA in dA_sequence:
        Cl_dA = roll_control_derivative(aileron_percent*bw/2,l1)
        Cl_p = roll_damping_coefficient(0,bw/2)
        P = steady_state_roll(Cl_dA,Cl_p,dA,V_stall)
        if P > 0:
            dt = time_req(45,P)
            if dt < 1.4:
                dA_upper = dA / (7/8)
                dA_lower = dA * (7/8)
                choices = pd.DataFrame({"ld":[aileron_lm_choice],
                                        "l1":[l1],
                                        "Ca_c":[Ca_c],
                                        "dA":[dA],
                                        "dA_upper":[dA_upper],
                                        "dA_lower":[dA_lower],
                                        "dt":[dt],
                                        "Steady roll rate":[P],
                                        "Cl_p":[Cl_p],
                                        "Cl_dA":[Cl_dA]})
                aileron_DP = aileron_DP.append(choices)
#aileron_DP = aileron_DP[(aileron_DP["dA_upper"] < 26.0) & (aileron_DP['dt'] < 1.39)].sort_values("dA_upper")
# print(aileron_DP)
aileron_final_design = aileron_DP.iloc[0]
aileron_lm = aileron_final_design.iloc[0]
aileron_l1 = aileron_final_design.iloc[1]
aileron_Ca_c = aileron_final_design.iloc[2]
aileron_dA = aileron_final_design.iloc[3]
aileron_dA_upper = aileron_final_design.iloc[4]
aileron_dA_lower = aileron_final_design.iloc[5]
aileron_dt = aileron_final_design.iloc[6]
aileron_P = aileron_final_design.iloc[7]
aileron_cl_p = aileron_final_design.iloc[8]
aileron_cl_dA = aileron_final_design.iloc[9]
# print(aileron_lm)
# print(aileron_l1)
# print(aileron_Ca_c)
# print(aileron_dA)
# print(aileron_dA_upper)
# print(aileron_dA_lower)
# print(aileron_dt)
# print(np.rad2deg(aileron_P))
M_cruise = V_cruise/a_cruise
c_accent_alpha_over_chalpha_theory_HT = 2.0/1.2 * 0.2
chalpha_theory = 6/2.5 * -0.2
Ch_alpha_bal_over_Chalpha = 0.85
c_accent_delta_over_chdelta_theory_HT = 0.8
chdelta_theory_HT = -0.72
Ch_delta_bal_over_Chdelta = 0.775
def Hingemoment_Coefficients_AF():
    c_accent_h_alpha= c_accent_alpha_over_chalpha_theory_HT*chalpha_theory
    Ch_alpha_bal = c_accent_h_alpha*Ch_alpha_bal_over_Chalpha
    Ch_alpha = Ch_alpha_bal/((1-M_cruise**2)**0.5)

    c_accent_h_delta = c_accent_delta_over_chdelta_theory_HT*chdelta_theory_HT
    ch_delta_bal = c_accent_h_delta*Ch_delta_bal_over_Chdelta
    Ch_delta = ch_delta_bal/((1-M_cruise**2)**0.5)
    #Ch_alpha = -0.02
    Ratio = Ch_alpha/Ch_delta
    return Ch_delta, Ch_alpha, Ratio , ch_delta_bal, Ch_alpha_bal

Eta_i = aileron_lm / (bw/2)
Eta_o = aileron_l1 / (bw/2)
Kalpha_i = 1.5
Kalpha_o = 3.6
DeltaCha_over_clalphaBK = 3.5 * 10**(-3)
B2 = 1
ch_alpha_M = Hingemoment_Coefficients_AF()[4] / ((1-M_cruise**2)**(1/2))
Sweep_hl = Calculate_wingsweep(0,0,0.76,taperw)
Kdelta_i = 1.4
Kdelta_o = 3.5
DeltaChd_cldBKd = 0.05
cl_delta = 4.8
ch_delta_M = Hingemoment_Coefficients_AF()[3] / ((1-M_cruise**2)**(1/2))
alpha_d = 0.52
def Hingemoment_Coefficients():
    Kalpha = Kalpha_i*(1-Eta_i)-Kalpha_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaC_h_alpha = DeltaCha_over_clalphaBK*(CL_alpha_AF*B2*Kalpha*np.cos(0))
    Ch_alpha = ((12*np.cos(0))/(12+2*np.cos(0)))*ch_alpha_M + DeltaC_h_alpha

    Kdelta = Kdelta_i*(1-Eta_i)-Kdelta_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaCh_delta =DeltaChd_cldBKd*(cl_delta*B2*Kdelta*np.cos(0)*np.cos(Sweep_hl))
    Ch_delta = np.cos(0)*np.cos(Sweep_hl)*(ch_delta_M +alpha_d*ch_alpha_M)*((2*np.cos(0))/(12+2*np.cos(0)))+DeltaCh_delta

    Ch_delta_horn = Ch_delta*0.26                                    # Horn Effect
    Ratio = Ch_alpha / Ch_delta_horn
    return Ch_delta_horn, Ch_alpha, Ratio, Ch_delta

l_top_lm = l_top(alpha_angle,beta_angle,aileron_lm,c_rw)
l_top_l1 = l_top(alpha_angle,beta_angle,aileron_l1,c_rw)
S_a = trapezoid_area(aileron_l1-aileron_lm,l_top_lm*0.3,l_top_l1*0.3)
y_m = (aileron_lm + aileron_l1)/2
#
control_F = -np.radians(aileron_dA)/0.4 * 1/2 * rho_5000 * V_stall * aileron_P * (bw/2) * S_a * aileron_Ca_c * ((Hingemoment_Coefficients()[1] * (2*y_m/bw)) - (Hingemoment_Coefficients()[3]*aileron_cl_p/(2*aileron_cl_dA)))
control_F_horn = -np.radians(aileron_dA)/0.4 * 1/2 * rho_5000 * V_stall * aileron_P * (bw/2) * S_a * aileron_Ca_c * ((Hingemoment_Coefficients()[1] * (2*y_m/bw)) - (Hingemoment_Coefficients()[0]*aileron_cl_p/(2*aileron_cl_dA)))
print(control_F)
print(control_F_horn)
# # # # print(aileron_cl_dA)
# # # # print(aileron_cl_p)
# # #print(-aileron_cl_dA/aileron_cl_p*np.radians(aileron_dA))
# # # print(aileron_l1)
# # # print(b/2)
# # # print(b/2 - 2.15/2)
