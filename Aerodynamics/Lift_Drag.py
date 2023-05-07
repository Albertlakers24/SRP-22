import numpy as np
import matplotlib.pyplot as plt
from Aerodynamics.AeroFunctions import *
from Aerodynamics.Drag import CD0_CR, CD0_TO, CD0_LD
from Constants.AircraftGeometry import bw, y_mac, taperw, Aw, l_f, DCL_Alpha0_TO, DCL_Alpha0_LD, delta_CL_TO, delta_CL_LD, CL_Alpha_TO, CL_Alpha_LD
from Constants.MissionInputs import M_cruise, Psi, phi, V_cruise, Calculate_Reynolds_number, kv_cruise

# PRELIMINARY
mach_TO = 0.3
beta_CR = Calculate_beta(M_cruise)
beta_TO = Calculate_beta(mach_TO)
sweep_half = Calculate_wingsweep(0, 0.25, 0.5, taperw)

''' AIRFOIL '''

Cl_list_AF = np.linspace(-1,1.3,100)
Alpha_list_AF = np.linspace(-12, 10, 100)
Cl_max_AF = 1.5
# Stall part
Alpha_list_stallAF = np.linspace(10,17, 40)
Cl_list_stallAF = []
for i in range(len(Alpha_list_stallAF)):
    x = Alpha_list_stallAF[i]
    Cl = -0.008*(x-15)**2 + 1.5
    Cl_list_stallAF.append(Cl)
Cl_list_stallAF = np.array(Cl_list_stallAF)

Cl_AF = np.append(Cl_list_AF, Cl_list_stallAF)
Alpha_AF = np.append(Alpha_list_AF, Alpha_list_stallAF)

''' LIFT - CRUISE '''
delta_alpha_CLMax = 2.1
Alpha0 = Alpha_list_AF[np.where(np.isclose(Cl_list_AF, 0, atol = 0.01))]#Alpha_AF_xfoil[np.where(np.isclose(CL_AF_xfoil,0, atol = 0.01))]

CL_max_CRw = Calculate_CL_max(Cl_max_AF,0)
CL_alpha_CRw = Calculate_CL_alpha(beta_CR,sweep_half)
Alpha_stall_CRw = Calculate_alpha_stall(CL_max_CRw, CL_alpha_CRw,np.radians(Alpha0), np.radians(delta_alpha_CLMax))
Alpha_trim = np.degrees(Calculate_alpha_trim(0.63, CL_alpha_CRw, np.radians(Alpha0)))

Alpha_list_CRw = np.arange(-10,10.25,0.25)
CL_list_CRw = []
for i in range(len(Alpha_list_CRw)):
    CL_value = CL_alpha_CRw * (np.radians(Alpha_list_CRw[i]) - np.radians(Alpha0))
    CL_list_CRw.append(CL_value)

Alpha_list_stallCRw = np.linspace(10,18, 40)
CL_list_stallCRw = []
for i in range(len(Alpha_list_stallCRw)):
    y = Alpha_list_stallCRw[i]
    #CL = -0.008*y**2 + 0.228*y - 0.186
    CL = -0.005*(y-15)**2 + 1.35
    CL_list_stallCRw.append(CL)

CL_list_stallCRw = np.array(CL_list_stallCRw)

CL_Wing_CR = np.append(CL_list_CRw, CL_list_stallCRw)
ALpha_Wing_CR = np.append(Alpha_list_CRw, Alpha_list_stallCRw)

'''DRAG - CRUISE '''

#- ----------------------- DRAG ---------------------------
h_winglet = 2.4
delta_A = 0#1.9*h_winglet/bw *Aw
A_eff = Aw + delta_A
twist_tip = -3
twist_root = 1
twist_list = np.linspace (twist_root, twist_tip, 50)
span_list = np.linspace(0, bw/2)
twist_mac = np.mean(twist_list[np.isclose(span_list, y_mac, rtol=0.05)]) #-0.74496644
delta_cd_twist = 0.00004*(twist_tip-twist_mac)
e_cruise = 1/(np.pi*Aw*Psi+(1/phi))

CD_list_CR = []
for i in range(len(CL_Wing_CR)):
    CDi_cruise = (CL_Wing_CR[i])**2/(np.pi *A_eff*e_cruise) + delta_cd_twist
    CD_cruise = (CD0_CR + CDi_cruise)
    CD_list_CR.append(CD_cruise)

CD_Wing_CR = np.array(CD_list_CR)

''' LIFT - TAKE OFF'''

CL0_TO = CL_Wing_CR[np.where(ALpha_Wing_CR==0)]
CL_list_TO = []
Alpha_Wing_TO = ALpha_Wing_CR + DCL_Alpha0_TO
Alpha0_TO = Alpha0 + DCL_Alpha0_TO
for i in range(len(Alpha_Wing_TO)-20):
    CL_TO = CL_Alpha_TO*(np.radians(Alpha_Wing_TO[i]) - np.radians(Alpha0_TO))
    CL_list_TO.append(CL_TO)
CL_list_TO = np.array(CL_list_TO)

#  STALL PART OF THE PLOT
Alpha_list_stallTO = np.linspace(8.1,14, 40)
CL_list_stall_TO = []
for i in range(len(Alpha_list_stallTO)):
    y = Alpha_list_stallTO[i]
    CL_TO = -0.0075*(y-12)**2 + 1.9
    CL_list_stall_TO.append(CL_TO)
CL_list_stall_TO = np.array(CL_list_stall_TO)

CL_list_TO= np.append(CL_list_TO, CL_list_stall_TO)
Alpha_Wing_TO = np.append(Alpha_Wing_TO[:-20], Alpha_list_stallTO)

'''DRAG - TAKE OFF'''
delta_f_TO = 15
Change_e = 0.0026 * delta_f_TO
e_TO = e_cruise + Change_e
CD_list_TO = []
for i in range(len(CL_list_TO)):
    CDi_TO = (CL_list_TO[i])**2/(np.pi *A_eff*e_TO) +delta_cd_twist
    CD_TO = (CD0_TO + CDi_TO)
    CD_list_TO.append(CD_TO)

CD_Wing_TO = np.array(CD_list_TO)


''' LIFT - LANDING '''

CL_list_LD = []
Alpha_Wing_LD = ALpha_Wing_CR + DCL_Alpha0_LD
Alpha0_LD = Alpha0 + DCL_Alpha0_LD
for i in range(len(Alpha_Wing_LD)-10):
    CL_LD = CL_Alpha_LD*(np.radians(Alpha_Wing_LD[i]) - np.radians(Alpha0_LD))
    CL_list_LD.append(CL_LD)
CL_list_LD = np.array(CL_list_LD)
#  STALL PART OF THE PLOT
Alpha_list_stallLD = np.linspace(7.15,11, 40)
CL_list_stall_LD = []
for i in range(len(Alpha_list_stallLD)):
    y = Alpha_list_stallLD[i]
    CL_LD = -0.0085*(y-9.5)**2 + 2.1
    CL_list_stall_LD.append(CL_LD)
CL_list_stall_LD = np.array(CL_list_stall_LD)

CL_list_LD= np.append(CL_list_LD, CL_list_stall_LD)
Alpha_Wing_LD = np.append(Alpha_Wing_LD[:-10], Alpha_list_stallLD)

''' DRAG - LANDING '''

delta_f_LD = 40
Change_e = 0.0026 * delta_f_LD
e_LD = e_cruise + Change_e
CD_list_LD = []
for i in range(len(CL_list_LD)):
    CDi_LD = (CL_list_LD[i])**2/(np.pi *A_eff*e_LD) + delta_cd_twist
    CD_LD = (CD0_LD + CDi_LD)
    CD_list_LD.append(CD_LD)

CD_Wing_LD = np.array(CD_list_LD)

''' AERODYNAMIC DATA ---------------------'''
CL_CD_CR = CL_Wing_CR/CD_Wing_CR
CL_CD_MaxCR = CL_CD_CR.max()
CL_CLCDMax = CL_Wing_CR[np.where(CL_CD_CR == CL_CD_MaxCR)]
Alpha_CLCDMax = ALpha_Wing_CR[np.where(CL_CD_CR == CL_CD_MaxCR)]

ind_cldes_CR = np.where(np.isclose(CL_Wing_CR, 0.63369832, atol = 0.01))
CL_Des_CR = CL_Wing_CR[ind_cldes_CR]        # CL Design Cruise
CD_Des_CR = CD_Wing_CR[ind_cldes_CR]        # CD Design Cruise
Alpha_Des_CR = ALpha_Wing_CR[ind_cldes_CR]         # Alpha Design Cruise
CL_CD_Des_CR = CL_CD_CR[ind_cldes_CR]

print('----------CRUISE DATA -----------')

print('CL/CD Cruise', CL_CD_Des_CR)
print('CL Cruise', CL_Des_CR)
print('CD Cruise', CD_Des_CR)
print('CL/CD Max Cruise', CL_CD_MaxCR)

''' TAKE-OFF ---------------------'''
CL_CD_TO = CL_list_TO/CD_Wing_TO
CL_CD_MaxTO = CL_CD_TO.max()

ind_cldes_TO = np.where(np.isclose(CL_list_TO,1.90, atol = 0.001))
CL_Des_TO = CL_list_TO[ind_cldes_TO]        # CL Design Cruise
CD_Des_TO = CD_Wing_TO[ind_cldes_TO]        # CD Design Cruise
Alpha_Des_TO = Alpha_Wing_TO[ind_cldes_TO]         # Alpha Design Cruise
CL_CD_Des_TO = CL_CD_TO[ind_cldes_TO]

print('----------TAKE OFF DATA -----------')

print('CL/CD Take-off', CL_CD_Des_TO)
print('CL Take-off', CL_Des_TO)
print('CD Take-off', CD_Des_TO)
print('CL/CD Max Tak-off', CL_CD_MaxTO)
print('Alpha Take off', Alpha_Des_TO)

''' LANDING ----------------------- '''
CL_CD_LD = CL_list_LD/CD_Wing_LD
CL_CD_MaxLD = CL_CD_LD.max()

ind_cldes_LD = np.where(np.isclose(CL_list_LD,2.11, atol = 0.01))
CL_Des_LD = CL_list_LD[ind_cldes_LD]        # CL Design Cruise
CD_Des_LD = CD_Wing_LD[ind_cldes_LD]        # CD Design Cruise
Alpha_Des_LD = Alpha_Wing_LD[ind_cldes_LD]         # Alpha Design Cruise
CL_CD_Des_LD = CL_CD_LD[ind_cldes_LD]
print('----------LANDING DATA -----------')
print('CL/CD Landing', CL_CD_Des_LD)
print('CL Landing', CL_Des_LD)
print('CD Landing', CD_Des_LD)
print('CL/CD Max Landing', CL_CD_MaxLD)
print('Alpha Landing', Alpha_Des_LD)

''' PLOTTING --------------------- '''
# plt.plot(Alpha_AF, Cl_AF, color = 'navy', linestyle = '-', label = 'Airfoil')
# plt.plot(ALpha_Wing_CR, CL_Wing_CR, color = 'coral', linestyle='-', label = 'Wing Clean')
# plt.plot(Alpha_Wing_TO, CL_list_TO, color = 'goldenrod', linestyle = '-', label = 'Wing @ Take-off')
# plt.plot(Alpha_Wing_LD, CL_list_LD, color = 'red', linestyle = '-', label = 'Wing @ Landing')
#
# plt.axhline(0,0,color = 'black')
# plt.axvline(0,0, color = 'black')
# plt.xlabel('Alpha[deg]')
# plt.ylabel('$C_l$[-]')
# plt.title('Cl-Alpha curve')
# plt.legend()
# plt.grid()
# plt.show()

