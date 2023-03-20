#from Constants.MissionInputs import *
import numpy as np
from Aerodynamics.AeroFunctions import Calculate_wingsweep
from Constants.AircraftGeometry import bw,taperw,c_rw, S_w, Aw, c_tw
from Constants.Aerodynamics import CL_Max_Clean, CL_MaxLand
def yield_func(low,high,step):
    i = low
    while i <= high:
        yield i
        i+=step
def HLD_TE_deltaClmax(Cf,df,flap_type):
    if flap_type == "single slotted":
        delta_c = -1/6600 * df *(df-104.2)      #Torenbeek page 533
        c_prime_over_c = 1 + delta_c*Cf
        delta_clmax_TE = 1.3
    if flap_type == "double slotted":
        if df < 15.01:
            delta_c = 0.3/15                                                    #Torenbeek page 533
        else:
            delta_c = 0.15 + df *((0.73-0.3)/(60-15))                   #Torenbeek page 533
        c_prime_over_c = 1 + delta_c * Cf
        delta_clmax_TE = 1.6 * c_prime_over_c
    if flap_type == "fowler":                                                   #Single slotted fowler
        if df < 10.01:
            delta_c = 0.45 / (1/1.2 * 10)                                       #Torenbeek page 533
        else:
            delta_c = 0.4 + df * (0.2/(45-(1/1.2*10)))            #Torenbeek page 533
        c_prime_over_c = 1 + delta_c*Cf
        delta_clmax_TE = 1.3 * c_prime_over_c

    return delta_clmax_TE,c_prime_over_c

def HLD_LE_deltaClmax(flap_type,c_prime_over_c):
    if flap_type == "LE_flap":
        delta_clmax_LE = 0.3
    if flap_type == "slat":
        delta_clmax_LE = 0.4 * c_prime_over_c
    else:
        delta_clmax_LE = 0
    return delta_clmax_LE

def wing_area(l_top,l_root,h):
    area = (l_top + l_root)/2 * h
    return area

def l_top(alpha,beta,height,base):
    x1 = height/ np.tan(alpha)
    x2 = height/ np.tan(beta)
    top = base - x1 - x2
    return top

#Design Options
te_hld = ["single slotted","double slotted","fowler"]
le_hld = ["None"]#["None","LE_flap","slat"]
lm = 1.5                                               #Minimum distance from center line to TE HLD [m]
lm_LE = 5.4                                            #Minimum distance from center line to LE HLD [m]

#Planform data
LE_hinge_line_angle_deg = Calculate_wingsweep(0,0.25,0.15,taperw)   #front spar position at 15% chord
TE_hinge_line_angle_deg = Calculate_wingsweep(0,0.25,0.60,taperw)    #rear spar position at 60% chord
LE_angle_deg = Calculate_wingsweep(0,0.25,0,taperw)                 #Leading edge angle
TE_angle_deg = Calculate_wingsweep(0,0.25,1,taperw)                 #Trailing edge angle
alpha_angle = np.radians(90 - LE_angle_deg)
beta_angle = np.radians(90 + TE_angle_deg)
aileron_percent = 0.55                                   #aileron starting position
# file = open("HLD_design_choices.txt", "w")
# file.write("ld[m]\tld[%]\tSwfS\tCf_C\tDCL\tDCLalpha\tDCl\tdf\tc_prime_over_c\tflap_type\n")
deltaCL_max_corrected = CL_MaxLand - CL_Max_Clean #Change it to CL_max_wing from airfoil selection

A1 = wing_area(c_tw,c_rw,bw/2)
C_t_trial = l_top(alpha_angle,beta_angle,bw/2,c_rw)
A1 = wing_area(C_t_trial,c_rw,bw/2)
print(C_t_trial)
print(c_tw)
print(A1)
print((c_rw+c_tw)/2 * bw/2)
print(S_w/2)
# print(l_top(alpha_angle,beta_angle,bw/2,c_rw))
# print(A1*2)
# print(LE_angle_deg)
# print(TE_angle_deg)
# print(A1)
# print(S_w/2)
# lm_top = l_top(alpha_angle,beta_angle,lm,c_rw)
# lm_topLE = l_top(alpha_angle,beta_angle,lm_LE,c_rw)
# l3_sequence = yield_func(lm_LE,(bw/2)*aileron_percent,0.01)
# for l3 in l3_sequence:
#     A2 = trapezoid_area(l3,l_top(alpha_angle,beta_angle,l3,c_rw),lm_topLE)
#     if 0.585 <= A2/(S_w/2):
#         ld = l3
#         ld_percent = l3/(bw/2) * 100
#         break
#     else:
#         ld = -1
#         ld_percent = -1
# print(lm_LE, " ", ld)
# for flap_type in te_hld:
#     for LE_type in le_hld:
#         deltaf = yield_func(10,40,5)
#         for df in deltaf:
#             SwfS_sequence = yield_func(0.3,0.9,0.01)
#             for SwfS in SwfS_sequence:
#                 Cf_sequnce = yield_func(0.10,0.41,0.01)
#                 for Cf in Cf_sequnce:
#                     DPS = False
#                     c_prime_over_c_LE_sequence = yield_func(1, 1.1, 0.01)
#                     for c_prime_over_c_LE in c_prime_over_c_LE_sequence:
#                         #Start of looped code
#                         deltaClmax_HLD, c_prime_over_c = HLD_TE_deltaClmax(Cf, df,flap_type)  # Calculating delta Cl max for trailing edge HLD
#                         deltaClmax_LE = HLD_LE_deltaClmax(LE_type,c_prime_over_c_LE)  # Calculating delta Cl max for leading edge HLD
#                         DP = 0.9 * deltaClmax_HLD * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg)) #+ 0.9 * deltaClmax_LE * 0.585 * np.cos(np.radians(LE_hinge_line_angle_deg))  # Calculating total CL max increase
#                         DPA = -15 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))  # Calculating change in CL alpha
#                         if DP >= deltaCL_max_corrected:
#                             l1_sequence = yield_func(lm,(bw/2)*aileron_percent,0.01)
#                             for l1 in l1_sequence:
#                                 A2 = trapezoid_area(l1, l_top(alpha_angle, beta_angle, l1, c_rw), lm_top)
#                                 if A2/ (S_w/2) >= SwfS
#                                     SwfS = A2/ (S_w/2)
#                                     DP = 0.9 * deltaClmax_HLD * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
#                                     DPA = -15 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
#                                     ld = l1
#                                     ld_percent = l1/(bw/2)*100
#                                     break
#                                 else:
#                                     ld = -1
#                                     ld_percent = -1
#                             file.write(str(round(ld,2)) + '\t' +str(round(ld_percent,2)) + '\t' +str(round(SwfS,2)) + '\t' + str(round(Cf,2))+ '\t' + str(round(DP,2))+ '\t' + str(round(DPA,2))+ '\t' + str(round(deltaClmax_HLD,2))+ '\t' + str(df)+ '\t' + str(round(c_prime_over_c,2))+ '\t' + flap_type + '+' + LE_type + '\n')
#                             DPS = True
#                         if DPS:    break
#                     if DPS:     break
# file.close()
# print("end code")