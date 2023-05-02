# from Aerodynamic_characteristics.AeroPlotsWing import CM_Wing, CM_AF_xfoil, CD_AF_xfoil, Alpha_AF_xfoil
# from Aerodynamic_characteristics.AeroPlotsTail import *
from Aerodynamics.Lift_Drag import ALpha_Wing_CR, Alpha_Wing_TO, Alpha_Wing_LD, CL_Wing_CR, CD_Wing_CR, CL_list_TO, CL_list_LD, Cl_AF, Alpha_AF, CD_list_TO, CD_list_LD, CD0_CR, CD0_LD, CD0_TO
import numpy as np
from Constants.Aerodynamics import CL_MaxLand
CD0_CR = CD0_CR
CD0_TO = CD0_TO
CD0_LD = CD0_LD

Alpha_Wing = ALpha_Wing
Alpha_Wing_TO = Alpha_Wing_TO
Alpha_Wing_LD = Alpha_Wing_LD

CL_list_cruise= CL_Wing_CR
CD_list_cruise= CD_Wing_CR
CL_list_TO= CL_list_TO
CD_list_TO= CD_list_TO
CL_list_LD= CL_list_LD
CD_list_LD = CD_list_LD

#CM_Wing = CM_Wing
CL_Wing = CL_Wing
CL_AF = Cl_AF
# CD_AF =CD_AF_xfoil
# CM_AF =CM_AF_xfoil
Alpha_AF = Alpha_AF
# Indices for slopes
ind_alpha1 =np.where(Alpha_Wing == -5)
ind_alpha2 =np.where(Alpha_Wing == 5)

# --------------------------- WING DATA ---------------------------

des_CL_CR = 0.6296         # SET - if approximation doesn't work - pls take CL from aero_data/Wing_634421.txt
des_CL_TO = CL_MaxLand/1.13         # delta CL = 0.62143 - delta f = 40 deg
des_CL_LD = CL_MaxLand #1.465         # delta CL = 0.74 - delta f = 40 deg

# print(np.where(np.isclose(CL_list_cruise, des_CL, rtol = 0.01)))
CL_CD_CR = CL_list_cruise/CD_Wing_CR
CL_CD_TO = CL_list_TO/CD_list_TO
CL_CD_LD = CL_list_LD/CD_list_LD

ind= np.where(np.isclose(CL_Wing_CR, 1, atol = 0.01))
CL_new = CL_Wing_CR[:71]
CD_new = CD_Wing_CR[: len(CL_new)]
#print(np.where(np.isclose(CL_datcom, 1, atol = 0.01)))
CL3_CD2 = CL_new**3/CD_new**2
print('max CL^3/CD^2', CL3_CD2.max())
ind_cl3_cd2 = np.where(CL3_CD2 == CL3_CD2.max())
print('CL at max CL^3/CD^2', CD_Wing_CR[ind_cl3_cd2])
print('CD at max CL^3/CD^2', CD_Wing_CR[ind_cl3_cd2])
# print(np.where(np.isclose(CL_list_cruise, des_CL_CR, rtol = 0.01)))
# Indices
def find_CL_CD(des_CL_CR):
    ind_clcdMax_CR = np.where(CL_CD_CR == CL_CD_CR.max())
    ind_cldes_CR = np.where(np.isclose(CL_list_cruise, des_CL_CR, atol = 0.01)) #np.where(CL_list_cruise == des_CL)
    CD_Des_CR = np.mean(CD_Wing_CR[ind_cldes_CR])
    return CD_Des_CR

ind_cldes_CR = np.where(np.isclose(CL_list_cruise, des_CL_CR, atol = 0.01)) #np.where(CL_list_cruise == des_CL)
CD_Des_CR = CD_Wing_CR[ind_cldes_CR]
# print(CD_Des_CR, find_CL_CD(des_CL_CR))

ind_clcdMax_CR = np.where(CL_CD_CR == CL_CD_CR.max())
ind_cldes_CR = np.where(np.isclose(CL_list_cruise, des_CL_CR, atol = 0.01)) #np.where(CL_list_cruise == des_CL)

ind_clcdMax_TO = np.where(CL_CD_TO == CL_CD_TO.max())
ind_cldes_TO = np.where(np.isclose(CL_list_TO, des_CL_TO, atol = 0.01)) #np.where(CL_list_cruise == des_CL)

ind_clcdMax_LD = np.where(CL_CD_LD == CL_CD_LD.max())
ind_cldes_LD = np.where(np.isclose(CL_list_TO, des_CL_LD, atol = 0.01)) #np.where(CL_list_cruise == des_CL)

ind_alpha0_Wing = np.where(Alpha_Wing == 0)

ind_cl_0_CR = np.where(np.isclose(CL_list_cruise, 0, atol = 0.01))
ind_cl_0_TO = np.where(np.isclose(CL_list_TO, 0, atol = 0.01))
ind_cl_0_LD = np.where(np.isclose(CL_list_TO, 0, atol = 0.01))

# USEFUL AERO DATA

# -------- CRUISE --------

CL_CLCDmax_CR = CL_list_cruise[ind_clcdMax_CR]  # CL @ CL/CD max Cruise
CD_CLCDmax_CR = CD_Wing_CR[ind_clcdMax_CR]  # CD @ CL/CD max Cruise
CL_CD_max_CR = CL_CD_CR[ind_clcdMax_CR]         # CL/CD max Cruise
Alpha_CLCDmax_CR = Alpha_Wing[ind_clcdMax_CR]   # Alpha @ CL/CD max Cruise

CL_Des_CR = CL_list_cruise[ind_cldes_CR]        # CL Design Cruise
CD_Des_CR = CD_Wing_CR[ind_cldes_CR]        # CD Design Cruise
Alpha_Des_CR = Alpha_Wing[ind_cldes_CR]         # Alpha Design Cruise
CL_CD_Des_CR = CL_CD_CR[ind_cldes_CR]           # CL/CD Design Cruise

CL0_wing_CR = CL_Wing[ind_alpha0_Wing]          # CL0 Wing - Cruise
CLmax_Wing_CR = CL_list_cruise.max()            # CL Max Wing - Cruise

AlphaCL0_CR = Alpha_Wing[ind_cl_0_CR]           # degrees

# print('CL @ CL/CD max Cruise -->', CL_CLCDmax_CR)
# print('CD @ CL/CD max Cruise -->', CD_CLCDmax_CR)
# print('CL/CD max Cruise -->', CL_CD_max_CR)
# # print('Alpha @ CL/CD max Cruise -->', Alpha_CLCDmax_CR)
# # #
print('CL Design Cruise -->', np.mean(CL_Des_CR))
print('CD @ CL Design Cruise -->', CD_Des_CR)
print('CL/CD @ CL Design Cruise -->', CL_CD_Des_CR)
print('Alpha @ CL Design Cruise -->', Alpha_Des_CR)
# #
# # print('CL_0 Wing Cruise -->', CL0_wing_CR)
print('CL_max Wing Cruise  -->', CLmax_Wing_CR)
# print('Alpha @ CL =0 Cruise', AlphaCL0_CR)
# -------- TAKE OFF --------

CL_CLCDmax_TO = CL_list_TO[ind_clcdMax_TO]      # CL @ CL/CD max Take-off
CD_CLCDmax_TO = CD_list_TO[ind_clcdMax_TO]      # CD @ CL/CD max Take-off
CL_CD_max_TO = CL_CD_TO[ind_clcdMax_TO]         # CL/CD max Take-off
Alpha_CLCDmax_TO = Alpha_Wing_TO[ind_clcdMax_TO]   # Alpha @ CL/CD max Take-off

CL_Des_TO = np.mean(CL_list_TO[ind_cldes_TO])            # CL Design Take-off
CD_Des_TO = CD_list_TO[ind_cldes_TO]            # CD Design Take-off
Alpha_Des_TO = Alpha_Wing_TO[ind_cldes_TO]         # Alpha Design Take-off
CL_CD_Des_TO = np.mean(CL_CD_TO[ind_cldes_TO] )          # CL/CD Design Take-off

CL0_wing_TO = CL_list_TO[ind_alpha0_Wing]       # CL0 Wing - Take-off
CLmax_Wing_TO = CL_list_TO.max()                # CL max Wing - Take-off

AlphaCL0_TO = Alpha_Wing_TO[ind_cl_0_TO]

# print('CL @ CL/CD max Take-off -->', CL_CLCDmax_TO)
# # print('CD @ CL/CD max Take-off -->', CD_CLCDmax_TO)
# print('CL/CD max Take-off -->', CL_CD_max_TO)
# # print('Alpha @ CL/CD max Take-off -->', Alpha_CLCDmax_TO)
# print('CL Design Take-off -->', CL_Des_TO)
# print('CD @ CL Design Take-off -->', CD_Des_TO)
# print('CL/CD @ CL Design Take-off -->', CL_CD_Des_TO)
# print('Alpha @ CL Design Take-off -->', Alpha_Des_TO)
# #
# print('CL_0 Wing Take-off -->', CL0_wing_TO)
# print('CL_max Wing Take-off -->', CLmax_Wing_TO)
# print('Alpha @ CL =0 Take-off', AlphaCL0_TO)

# -------- Landing --------

CL_CLCDmax_LD = CL_list_TO[ind_clcdMax_LD]      # CL @ CL/CD max Landing
CD_CLCDmax_LD = CD_list_LD[ind_clcdMax_LD]      # CD @ CL/CD max Landing
CL_CD_max_LD = CL_CD_LD[ind_clcdMax_LD]         # CL/CD max Landing
Alpha_CLCDmax_LD = Alpha_Wing[ind_clcdMax_LD]   # Alpha @ CL/CD max Landing

CL_Des_LD = CL_list_TO[ind_cldes_LD]            # CL Design Landing
CD_Des_LD = CD_list_LD[ind_cldes_LD]            # CD Design Landing
Alpha_Des_LD = Alpha_Wing[ind_cldes_LD]         # Alpha Design Landing
CL_CD_Des_LD = CL_CD_LD[ind_cldes_LD]           # CL/CD Design Landing

CL0_wing_LD = CL_list_TO[ind_alpha0_Wing]       # CL0 Wing - Landing
CLmax_Wing_LD = CL_list_TO.max()                # CL max Wing - Landing
AlphaCL0_LD = Alpha_Wing[ind_cl_0_LD]

ind_alpha_Wing1_TO= np.where(Alpha_Wing == -5)
ind_alpha_Wing2_TO = np.where(Alpha_Wing == 5)

CL_alpha_TO =((CL_list_TO[ind_alpha_Wing2_TO] - CL_list_TO[ind_alpha_Wing1_TO])/10)*180/np.pi    # CL_alpha Wing slope
# print('Cm_alpha Wing -->', Cm_alpha)
# print('CL_alpha Wing at Take-off -->', CL_alpha_TO)

# print('CL @ CL/CD max Landing -->', CL_CLCDmax_LD)
# # print('CD @ CL/CD max Landing -->', CD_CLCDmax_LD)
# print('CL/CD max Landing -->', CL_CD_max_LD)
# print('Alpha @ CL/CD max Landing -->', Alpha_CLCDmax_LD)
#
# print('CL Design Landing -->', CL_Des_LD)
# print('CD @ CL Design Landing -->', CD_Des_LD)
# print('CL/CD @ CL Design Landing -->', CL_CD_Des_LD)
# print('Alpha @ CL Design Landing -->', Alpha_Des_LD)

# print('CL_0 Wing Landing -->', CL0_wing_LD)
# print('CL_max Wing Landing -->', CLmax_Wing_LD)
# print('Alpha @ CL =0 Landing', AlphaCL0_LD)

# -------- Common --------


#Cm_alpha =(CM_Wing[ind_alpha2] - CM_Wing[ind_alpha1])/10                # [deg^-1] Cm_alpha Wing slope
CL_alpha =((CL_Wing[ind_alpha2] - CL_Wing[ind_alpha1])/10)*180/np.pi    # [rad^-1] CL_alpha Wing slope
# print('Cm_alpha Wing -->', Cm_alpha)
# print('CL_alpha Wing -->', CL_alpha)

# delta_alphaCL0_TO = AlphaCL0_CR - AlphaCL0_TO
# delta_alphaCL0_LD = AlphaCL0_CR - AlphaCL0_LD

# print('Delta Alpha@CL0 (Cruise - Take-off', delta_alphaCL0_TO)
# print('Delta Alpha@CL0 (Cruise - Landing', delta_alphaCL0_LD)

# ---------------------- WING airfoil ----------------------

#CL_CD_AF = CL_AF/CD_AF

# Indices
ind_alpha0_AF = np.where(np.isclose(Alpha_AF, 0, atol = 0.1))
ind_alpha_AFWing1 = np.where(np.isclose(Alpha_AF, -5, atol = 1))
ind_alpha_AFWing2 = np.where(np.isclose(Alpha_AF, 5, atol = 1))

# USEFUL DATA
#Cm0_AF = CM_AF[np.where(np.isclose(Alpha_AF, 0, atol = 0.01))]                     # Cm0 airfoil
Cl_alpha_AF =np.mean(((CL_AF[ind_alpha_AFWing2] - CL_AF[ind_alpha_AFWing1])/10)) *180/np.pi # Cl_alpha wing airfoil slope

# print('Cm_0 - wing airfoil -->', Cm0_AF)
# print('Cl_alpha - wing airfoil -->', Cl_alpha_AF)
# ind_CD0_CLDes = np.where(np.isclose(Alpha_AF_xfoil, Alpha_Des_TO.mean(), rtol = 0.01))
# Cd0_AF_CLDes = CD_AF[ind_CD0_CLDes]       # CD0 airfoil at CL Design
# Cd0_AF_Max = CD_AF.max()
# Cd0_AF_ALpha0 = CD_AF[np.where(Alpha_AF_xfoil == 0)]
# # print('Cd0 airfoil @ CL des', Cd0_AF_CLDes)
# print('CD0 @ Alpha = 0', Cd0_AF_ALpha0)
# print(CD_AF_xfoil.min())

# ---------------------- TAIL ----------------------

# ind_alpha_AFTail1 = np.where(Alpha_AF_Tail == -5)
# ind_alpha_AFTail2 = np.where(Alpha_AF_Tail == 5)
# ind_alphaTailCruise = np.where(Alpha_AF_Tail == 0.5)

# # USEFUL DATA for AIRFOIL (0012)
# Cl_alpha_AFtail =((CL_AF_Tail[ind_alpha_AFTail2] - CL_AF_Tail[ind_alpha_AFTail1])/10)*180/np.pi # Cl_alpha tail airfoil slope
#
# # USEFUL DATA for Horizontal Tail
#
# ind_alpha_HTail1 = np.where(Alpha_HTail == -5)
# ind_alpha_HTail2 = np.where(Alpha_HTail == 5)
#
# CL_alpha_htail =((CL_HTail[ind_alpha_HTail2] - CL_HTail[ind_alpha_HTail1])/10)*180/np.pi    # [rad^-1] Cl_alpha Horizontal tail slope
# Cm_alpha_htail =(CM_HTail[ind_alpha_HTail2] - CM_HTail[ind_alpha_HTail1])/10                # Cm_alpha Wing slope
# CL_HTail_CR = CL_HTail[ind_alphaTailCruise]
# # print(Cm_alpha_htail, Cm_alpha)
# # USEFUL DATA for Vertical Tail
# # print('CL H Tail Cruise', CL_HTail_CR)
# ind_alpha_VTail1 = np.where(Alpha_VTail == -5)
# ind_alpha_VTail2 = np.where(Alpha_VTail == 5)
#
# CL_alpha_vtail =((CL_VTail[ind_alpha_VTail2] - CL_VTail[ind_alpha_VTail1])/10)*180/np.pi   # Cl_alpha Vertical tail slope
# #
# print('Cl_alpha - tail airfoil -->', Cl_alpha_AFtail)
# # print('CL_alpha - H tail -->', CL_alpha_htail)
# # print('CL_alpha - V tail -->', CL_alpha_vtail)
#
# print('Cl alpha - tail airfoil', 0.1*180/np.pi)


