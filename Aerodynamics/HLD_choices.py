#from Constants.MissionInputs import *
import numpy as np
from Aerodynamics.AeroFunctions import Calculate_wingsweep
from Constants.AircraftGeometry import bw,taperw,c_rw, S_w, Aw, c_tw
from Constants.Aerodynamics import CL_Max_Clean, CL_MaxLand, CL_Alpha_Wing
import pandas as pd
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
            delta_c = 0.45 / (1/1.2 * 10) * df                                      #Torenbeek page 533
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

def trapezoid_area(h,top,base):
    area = (top+base)/2 * h
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
TE_hinge_line_angle_deg = Calculate_wingsweep(0,0.25,0.55,taperw)    #rear spar position at 60% chord
LE_angle_deg = Calculate_wingsweep(0,0.25,0,taperw)                 #Leading edge angle
TE_angle_deg = Calculate_wingsweep(0,0.25,1,taperw)                 #Trailing edge angle
alpha_angle = np.radians(90 - LE_angle_deg)
beta_angle = np.radians(90 + TE_angle_deg)
aileron_percent = 0.55                                   #aileron starting position
data = []
deltaCL_max_corrected = CL_MaxLand - CL_Max_Clean #Change it to CL_max_wing from airfoil selection
lm_top = l_top(alpha_angle,beta_angle,lm,c_rw)
for flap_type in te_hld:
    for LE_type in le_hld:
        deltaf = yield_func(10,45,1)
        for df in deltaf:
            SwfS_sequence = yield_func(0.3,0.9,0.01)
            for SwfS in SwfS_sequence:
                Cf_sequnce = yield_func(0.10,0.41,0.01)
                for Cf in Cf_sequnce:
                    DPS = False
                    c_prime_over_c_LE_sequence = yield_func(1, 1.1, 0.01)
                    for c_prime_over_c_LE in c_prime_over_c_LE_sequence:
                        #Start of looped code
                        deltaClmax_HLD, c_prime_over_c = HLD_TE_deltaClmax(Cf, df,flap_type)  # Calculating delta Cl max for trailing edge HLD
                        deltaClmax_LE = HLD_LE_deltaClmax(LE_type,c_prime_over_c_LE)  # Calculating delta Cl max for leading edge HLD
                        DP = 0.9 * deltaClmax_HLD * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg)) #+ 0.9 * deltaClmax_LE * 0.585 * np.cos(np.radians(LE_hinge_line_angle_deg))  # Calculating total CL max increase
                        DPA = -15 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))  # Calculating change in CL alpha
                        if DP >= deltaCL_max_corrected:
                            l1_sequence = yield_func(lm,(bw/2)*aileron_percent,0.01)
                            for l1 in l1_sequence:
                                A2 = trapezoid_area(l1, l_top(alpha_angle, beta_angle, l1, c_rw), lm_top)
                                if A2/ (S_w/2) >= SwfS:
                                    SwfS = A2/ (S_w/2)
                                    DP = 0.9 * deltaClmax_HLD * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
                                    DPA = -15 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
                                    ld = l1
                                    ld_percent = l1/(bw/2)*100
                                    break
                                else:
                                    ld = -1
                                    ld_percent = -1
                            data.append([ld,ld_percent,SwfS,Cf,DP,DPA,deltaClmax_HLD,df,c_prime_over_c,flap_type])
                            DPS = True
                        if DPS:    break
                    if DPS:     break
HLD_choices = pd.DataFrame(data)
HLD_choices.columns = ["ld[m]","ld[%]","SwfS","Cf_C","DCL","DCLalpha","DCl","df","c_prime_over_c","flap_type"]
HLD_choices_single = HLD_choices[(HLD_choices['ld[m]']>0) & (HLD_choices["flap_type"] == "single slotted")].sort_values("ld[%]")
HLD_choices_fowler = HLD_choices[(HLD_choices['ld[m]']>0) & (HLD_choices["flap_type"] == "fowler") & (HLD_choices['df']  == 40) & (HLD_choices['Cf_C'] < 0.37) & (HLD_choices['Cf_C'] > 0.33)].sort_values("ld[%]")
HLD_final_design = HLD_choices_fowler.iloc[0]
ld = HLD_final_design[0]
ld_perc = HLD_final_design[1]
SwfS = HLD_final_design[2]
Cf_C = HLD_final_design[3]
DCL = HLD_final_design[4]
DPA = HLD_final_design[5]
DCl = HLD_final_design[6]
flap_deflection = HLD_final_design[7]
c_prime_over_c_land = HLD_final_design[8]
S_total_land = S_w * (1 + SwfS*(c_prime_over_c_land-1))
cf_over_cprime = Cf_C * (1/c_prime_over_c_land)
c_prime_over_c_to = (0.4 + 15 * (0.2/(45-(1/1.2*10))))*Cf_C + 1
S_total_to = S_w * (1 + SwfS*(c_prime_over_c_to-1))
DCL_to = 0.9 * HLD_TE_deltaClmax(Cf_C,20,"fowler")[0] * 0.6 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
DPA_to = -10 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
DCL_alpha_land = S_total_land/S_w * CL_Alpha_Wing
DCL_alpha_to = S_total_to/S_w * CL_Alpha_Wing