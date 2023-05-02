import numpy as np
import matplotlib.pyplot as plt
from Aerodynamics.AeroFunctions import *
from Aerodynamics.Drag import CD0_CR, CD0_TO, CD0_LD
from Constants.AircraftGeometry import bw, y_mac, taperw, Aw, l_f
from Constants.MissionInputs import M_cruise, Psi, phi, V_cruise, Calculate_Reynolds_number, kv_cruise

mach_TO = 0.3
# -------------- CL - alpha curve approx ---------------------
Cl_list_theory = np.linspace(-1,1.3,100)
Alpha_list_theoryAF = np.linspace(-12, 10, 100)

#Cl_list_theory_stall = np.append(np.linspace(1.3,1.5, 20), np.linspace(1.5, 1.3, 20))
Alpha_list_theory_stallAF = np.linspace(10,17, 40)
Cl_list_theory_stall = []
for i in range(len(Alpha_list_theory_stallAF)):
    x = Alpha_list_theory_stallAF[i]
    Cl = -0.008*(x-15)**2 + 1.5
    Cl_list_theory_stall.append(Cl)

Cl_list_theory_stall = np.array(Cl_list_theory_stall)
# ---------------------- DATCOM METHOD ----------------------

delta_alpha_CLMax = 2.1 #25.46 * 0.21
Alpha0 = Alpha_list_theoryAF[np.where(np.isclose(Cl_list_theory, 0, atol = 0.01))]#Alpha_AF_xfoil[np.where(np.isclose(CL_AF_xfoil,0, atol = 0.01))]

Cl_max_AF = 1.5 #CL_AF_xfoil.max()
beta_CR = Calculate_beta(M_cruise)
beta_TO = Calculate_beta(mach_TO)
sweep_half = Calculate_wingsweep(0, 0.25, 0.5, taperw)

CL_max_datcom = Calculate_CL_max(Cl_max_AF,0)
CL_alpha_datcom = Calculate_CL_alpha(beta_CR,sweep_half)
CL_alpha_TO_clean = Calculate_CL_alpha(beta_TO, sweep_half)
Alpha_stall_datcom = Calculate_alpha_stall(CL_max_datcom, CL_alpha_TO_clean,np.radians(Alpha0), np.radians(delta_alpha_CLMax))
Alpha_trim = np.degrees(Calculate_alpha_trim(0.63, CL_alpha_datcom, np.radians(Alpha0)))
Alpha_datcom = np.arange(-10,10.25,0.25)
CL_datcom = np.array([])
for i in range(len(Alpha_datcom)):
    CL_value = CL_alpha_datcom * (np.radians(Alpha_datcom[i]) - np.radians(Alpha0))
    CL_datcom = np.append(CL_datcom, CL_value)
print(Alpha0)
print('CL alpha', CL_alpha_datcom)
print('Alpha stall', Alpha_stall_datcom, 'rad',  np.degrees(Alpha_stall_datcom), 'deg')
print('CL max', CL_max_datcom)
print('Cl_datcom final', CL_datcom[-1], Alpha_datcom[-1])
print('print value', Alpha0)
print('delta Cl max alpha', delta_alpha_CLMax)
print('Alpha trim', Alpha_trim)
print('Alpha 0', Alpha0)

Alpha_list_theory_stall = np.linspace(10,18, 40)
CL_list_theory_stall = []
for i in range(len(Alpha_list_theory_stall)):
    y = Alpha_list_theory_stall[i]
    #CL = -0.008*y**2 + 0.228*y - 0.186
    CL = -0.005*(y-15)**2 + 1.35
    CL_list_theory_stall.append(CL)

CL_list_theory_stall = np.array(CL_list_theory_stall)
print(Alpha_list_theory_stall[0], Alpha_datcom[-1])

CL_Wing = np.append(CL_datcom, CL_list_theory_stall)
ALpha_Wing = np.append(Alpha_datcom, Alpha_list_theory_stall)

Cl_AF = np.append(Cl_list_theory, Cl_list_theory_stall)
Alpha_AF = np.append(Alpha_list_theoryAF, Alpha_list_theory_stallAF)

# ------------------ HLDs ------------------

#DCL = 0.94
DCL_Alpha0_TO =-6.2
DCL_Alpha0_LD =-9.3
red_CD = 0.9
dCL_TO = 0.9
dCL_LD = 1.05
dCL_alpha_TO = 1.1
dCL_alpha_LD = 1.1495
ClAlpha = CL_alpha_datcom
CL_Alpha_TO = 6.14#dCL_alpha_TO * ClAlpha
CL_Alpha_LD = 6.20 #dCL_alpha_LD * ClAlpha

Alpha_0 = ALpha_Wing[np.where(ALpha_Wing==0)]
Alpha_0_TO = Alpha_0 + DCL_Alpha0_TO
Alpha_0_LD = Alpha_0 + DCL_Alpha0_LD


#CL_list_cruise = CL_Wing
# CL_list_TO = CL_Wing + dCL_TO
# CL_list_LD = CL_Wing + dCL_LD
# ------------- TAKE OFF ----------------------------------
CL_list_TO = [] #CL_Wing + dCL_TO

Alpha_Wing_TO = ALpha_Wing + DCL_Alpha0_TO

for i in range(len(Alpha_Wing_TO)-10):
    dCL = CL_Alpha_TO*(np.radians(ALpha_Wing[i]) - np.radians(Alpha_0)) + CL_Wing[np.where(ALpha_Wing==0)]
    CL_list_TO.append(dCL)
CL_list_TO = np.array(CL_list_TO)

print('CL alpha TO', CL_Alpha_TO)
print(CL_list_TO[-1], Alpha_Wing_TO[-10])
#CL_list_LD = CL_Wing + dCL_LD

Alpha_stall_TO = Calculate_alpha_stall(2.25, CL_Alpha_TO,np.radians(Alpha_0_TO), np.radians(delta_alpha_CLMax))
print('Alpha stall TO', np.degrees(Alpha_stall_TO))

#  STALL PART OF THE PLOT
Alpha_list_stall_TO = np.linspace(9.4,16, 40)
CL_list_stall_TO = []
for i in range(len(Alpha_list_stall_TO)):
    y = Alpha_list_stall_TO[i]
    #CL = -0.008*y**2 + 0.228*y - 0.186
    CL_TO = -0.018*(y-13.8)**2 + 2.25
    CL_list_stall_TO.append(CL_TO)
CL_list_stall_TO = np.array(CL_list_stall_TO)

CL_list_TO= np.append(CL_list_TO, CL_list_stall_TO)
Alpha_Wing_TO = np.append(Alpha_Wing_TO[:-10], Alpha_list_stall_TO)
print("CL max",CL_list_TO.max())
# ---------------------- LANDING - FLAPS 40 DEG --------------
CL_list_LD = [] #CL_Wing + dCL_TO

Alpha_Wing_LD = ALpha_Wing + DCL_Alpha0_LD

for i in range(len(Alpha_Wing_LD)):
    dCL = CL_Alpha_LD*(np.radians(ALpha_Wing[i]) - np.radians(Alpha_0)) + CL_Wing[np.where(ALpha_Wing==0)]
    CL_list_LD.append(dCL)
CL_list_LD = np.array(CL_list_LD)

print('CL alpha landing', CL_Alpha_LD)
print(CL_list_LD[-1], Alpha_Wing_LD[-1])

#CL_list_LD = CL_Wing + dCL_LD

Alpha_stall_LD = Calculate_alpha_stall(CL_max_datcom, CL_Alpha_LD ,np.radians(Alpha_0_LD), np.radians(delta_alpha_CLMax))
print('alpha stall landing', np.degrees(Alpha_stall_LD))

#  STALL PART OF THE PLOT
Alpha_list_stall_LD = np.linspace(8.7,16, 40)
CL_list_stall_LD = []
for i in range(len(Alpha_list_stall_LD)):
    y = Alpha_list_stall_LD[i]
    #CL = -0.008*y**2 + 0.228*y - 0.186
    CL_LD = -0.015*(y-12.5)**2 + 2.4
    CL_list_stall_LD.append(CL_LD)
CL_list_stall_TO = np.array(CL_list_stall_TO)

CL_list_LD= np.append(CL_list_LD, CL_list_stall_LD)
Alpha_Wing_LD = np.append(Alpha_Wing_LD, Alpha_list_stall_LD)
print("CL max",CL_list_LD.max())



#- ----------------------- DRAG ---------------------------
A = Aw     #Aspect Ratio (12-14) #Reference to ATR 72
h_winglet = 2.4
delta_A = 1.9*h_winglet/bw *A
A_eff = A + delta_A
twist_tip = -3
twist_mac =  -0.74496644
delta_cd_twist = 0.00004*(twist_tip-twist_mac)
e_cruise = 1/(np.pi*A*Psi+(1/phi))
print('e', e_cruise)
CD0_Cruise = CD0_CR
CD_list_cruise = []
for i in range(len(CL_Wing)):
    CDi_cruise = (CL_Wing[i])**2/(np.pi *A_eff*e_cruise) + delta_cd_twist
    CD_cruise = (CD0_Cruise + CDi_cruise)
    CD_list_cruise.append(CD_cruise)


delta_f_TO = 15
Change_e = 0.0026 * delta_f_TO
e_TO = e_cruise + Change_e
CD_list_TO = []
for i in range(len(CL_list_TO)):
    CDi_TO = (CL_list_TO[i])**2/(np.pi *A_eff*e_TO) +delta_cd_twist
    CD_TO = (CD0_TO + CDi_TO)
    CD_list_TO.append(CD_TO)

delta_f_LD = 40
Change_e = 0.0026 * delta_f_LD
e_LD = e_cruise + Change_e
CD_list_LD = []
for i in range(len(CL_list_LD)):
    CDi_LD = (CL_list_LD[i])**2/(np.pi *A_eff*e_LD) + delta_cd_twist
    CD_LD = (CD0_LD + CDi_LD)
    CD_list_LD.append(CD_LD)

CD_list_cruise = np.array(CD_list_cruise)
CD_list_TO = np.array(CD_list_TO)
CD_list_LD = np.array(CD_list_LD)


# ----------------- TAIL AIRFOIL --------------------------
# -------------- CL - alpha curve approx ---------------------
Cl_list_tail = np.linspace(-1.57,1.57,100)
Alpha_list_tailAF = np.linspace(-15.4, 15.4, 100)
print(Cl_list_tail[-1], Alpha_list_tailAF[-1])
#Cl_list_theory_stall = np.append(np.linspace(1.3,1.5, 20), np.linspace(1.5, 1.3, 20))
Alpha_list_tailstall_AF = np.linspace(15.5,17, 40)
Cl_list_tail_stall = []
for i in range(len(Alpha_list_tailstall_AF)):
    x = Alpha_list_tailstall_AF[i]
    #Cl = -0.02*x**2 + 0.54*x - 2.1
    Cl = -0.1*(x-16)**2 + 1.6
    Cl_list_tail_stall.append(Cl)

print((Cl_list_tail[20] - Cl_list_tail[10])/(Alpha_list_tailAF[20] - Alpha_list_tailAF[10]) * 180/np.pi)
print((Cl_list_tail[40] - Cl_list_tail[30])/(Alpha_list_tailAF[40] - Alpha_list_tailAF[30]))
print(Calculate_Reynolds_number(V_cruise, l_f, kv_cruise))

Alpha_list_tailAF_neg = np.linspace(-17, -15.5, 40)
Cl_list_tail_stall_neg = []
for i in range(len(Alpha_list_tailAF_neg)):
    x = Alpha_list_tailAF_neg[i]
    #Cl = -0.02*x**2 + 0.54*x - 2.1
    Cl = 0.1*(x+16)**2 - 1.6
    Cl_list_tail_stall_neg.append(Cl)

Cl_list_tail_stall = np.array(Cl_list_tail_stall)
Cl_list_tail_stall_neg = np.array(Cl_list_tail_stall_neg)
Cl_AF_tail = np.append(Cl_list_tail_stall_neg, Cl_list_tail)
Cl_AF_tail = np.append(Cl_AF_tail, Cl_list_tail_stall)
Alpha_AF_tail = np.append(Alpha_list_tailAF_neg, Alpha_list_tailAF)
Alpha_AF_tail = np.append(Alpha_AF_tail, Alpha_list_tailstall_AF)
# ------------------ CL-ALPHA PLOTTING ------------------

# plt.plot(Alpha_AF, Cl_AF, color = 'navy', linestyle = '-', label = 'Airfoil')
# plt.plot(ALpha_Wing, CL_Wing, color = 'coral', linestyle='-', label = 'Wing')
# plt.plot(Alpha_Wing_TO, CL_list_TO, color = 'goldenrod', linestyle = '-', label = 'Wing @ Flaps = 15 deg')
# plt.plot(Alpha_Wing_LD, CL_list_LD, color = 'teal', linestyle = '-', label = 'Wing @ Flaps = 40 deg')
# plt.axhline(0,0,color = 'black')
# plt.axvline(0,0, color = 'black')
# plt.xlabel('Alpha[deg]')
# plt.ylabel('$C_l$[-]')
# plt.title('Cl-Alpha curve')
# plt.legend()
# plt.grid()
# plt.show()

#
# # ----------------- CL-CD PLOTS ----------------
# plt.plot(CD_list_cruise, CL_Wing , color = 'coral', linestyle='-', label = 'Wing - clean')
# plt.plot(CD_list_TO, CL_list_TO, color = 'goldenrod', linestyle = '-', label = 'Wing @ Flaps = 15 deg')
# plt.plot(CD_list_LD, CL_list_LD, color = 'teal', linestyle = '-', label = 'Wing @ Flaps = 40 deg')
# plt.axhline(0,0,color = 'black')
# plt.axvline(0,0, color = 'black')
# plt.xlabel('$C_D$[-]')
# plt.ylabel('$C_L$[-]')
# plt.title('Drag Polar')
# plt.legend()
# plt.grid()
# plt.show()

# # ------------------ CD-ALPHA PLOTTING ------------------
#
# #plt.plot(Alpha_AF, Cl_AF, color = 'navy', linestyle = '--', label = 'Airfoil')
# plt.plot(ALpha_Wing, CD_list_cruise, color = 'coral', linestyle='-', label = 'Wing')
# plt.plot(Alpha_Wing_TO, CD_list_TO, color = 'goldenrod', linestyle = '-', label = 'Wing @ Flaps = 15 deg')
# plt.plot(Alpha_Wing_LD, CD_list_LD, color = 'teal', linestyle = '-', label = 'Wing @ Flaps = 40 deg')
#
# plt.xlabel('Alpha[deg]')
# plt.ylabel('CD[-]')
# plt.title('CD-Alpha curve')
# plt.legend()
# plt.grid()
# plt.show()


# ----------- TAIL plotting
# plt.plot(Alpha_AF_tail, Cl_AF_tail, color = 'mediumvioletred')
# plt.axhline(0,0,color = 'black')
# plt.axvline(0,0, color = 'black')
# plt.xlabel('Alpha[deg]')
# plt.ylabel('Cl[-]')
# plt.title('Cl-Alpha curve - Tail airfoil')
# plt.grid()
# plt.show()

a = np.linspace(0, bw/2, 150)
c = np.linspace(3.2171476 ,1.4477164,150)
ind = np.isclose(a, 3.01/2, rtol=0.05)
print(a[ind])
print(c[ind])
