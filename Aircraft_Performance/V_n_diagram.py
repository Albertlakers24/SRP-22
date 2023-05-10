import numpy as np
import matplotlib.pyplot as plt
from Constants.MissionInputs import rho_0,rho_cruise,lbs_kg,V_cruise, ISA_calculator, h_cruise, dt_cruise,dt_takeoff,h_loiter,dt_loiter,g,Aw
from Constants.Aerodynamics import CL_Max_Clean,CL_Alpha_Wing
from Constants.Masses_Locations import m_mto,W_S_design, m_oem, m_zf, beta_s_land_fc
from Constants.AircraftGeometry import S_w, bw, c_mac_w

# Density
rho = 0 # just to define rho, this NEVER changes
density = 0 # 0 @ sea-level, 1 @ cruise, else @ loiter
if density == 0:
    rho = rho_0
elif density == 1:
    rho = rho_cruise
else:
    rho = 0.663838

if density == 0:
    rho = rho_0 # density at sea level
elif density == 1:
    rho = rho_cruise
else:
    rho = 0.663838

#----------------------------------------MANEUVER DIAGRAM DESIGN-------------------------------------------------------
#Max lift coefficient
if density == 0:
    Clmax = 1.45 #CL_max_landing
elif density == 1:
    Clmax = CL_Max_Clean
else:
    Clmax = CL_Max_Clean # add loiter Clmax to Constants file

# print(Clmax)

#Load factor
m_mto_lbs = m_mto / lbs_kg #conversion between kg and lbs
n = 2.1 + (24000 / (m_mto_lbs + 10000))

if n < 2.5:
    n_max = 2.5
elif n > 3.8:
    n_max = 3.8
else:
    n_max = n
# print(n_max)
n_min = -1

#Cruise speed EAS
V_C = V_cruise
V_C = V_C * np.sqrt(rho/rho_0)

#Stall speed EAS
V_S = np.sqrt((2 * 1 * W_S_design) / (rho * Clmax)) * np.sqrt(rho / rho_0)

#Dive speed EAS
if density == 0:
    h = 0
    dt = dt_takeoff
elif density == 1:
    h = h_cruise
    dt = dt_cruise
else:
    h = h_loiter
    dt = dt_loiter

a = ISA_calculator(h,dt)[3]

isa_cruise = ISA_calculator(h_cruise, dt_cruise)
# print("ISA CRUISE", isa_cruise)

M_C = V_C / a
M_D = M_C + 0.05

V_D1 = V_C/0.8
V_D2 = M_D * a
V_D = min(V_D1, V_D2)
#---- CONSTRUCTION OF MANEUVER DIAGRAM ----
#From 0 to V_A
V_A = []
n_A_list = []
V_A_fin = V_S * np.sqrt(n_max)
for i in np.arange(0, V_A_fin, 0.1):
    V_A.append(i)
    n = 0.5 * i**2 * rho_0 * Clmax / W_S_design
    n_A_list.append(n)

#From 0 to V_H
V_H = []
n_H_list = []
for j in np.arange(0, V_S, 0.1):
    V_H.append(j)
    n = -0.5 * j**2 * rho_0 * Clmax / W_S_design
    n_H_list.append(n)

# FLAPS
V_max_flaps = np.sqrt(2 * W_S_design / (0.5 * rho_0 * 2.4))
V_intersect = 79.6987
V_A_flaps = []
n_A_list_flaps = []
V_A_fin_flaps = V_S * np.sqrt(n_max)
for i in np.arange(0, V_max_flaps, 0.1):
    V_A_flaps.append(i)
    n = 0.5 * i**2 * rho_0 * 2.4 / W_S_design
    n_A_list_flaps.append(n)

# print(n_A_list_flaps)

#Graph
point_F = [V_C, -1]
point_E = [V_D, 0]
point_A = [V_A_fin, n_max]
point_H = [V_S, -1]
x_values_FE = [point_F[0], point_E[0]]
y_values_FE = [point_F[1], point_E[1]]
point_D = [V_D, n_max]
x_values_DE = [point_D[0], point_E[0]]
y_values_DE = [point_D[1], point_E[1]]
x_values_AD = [point_A[0], point_D[0]]
y_values_AD = [point_A[1], point_D[1]]

# fig, ax = plt.subplots()
# plt.tick_params(labelbottom = False)
# plt.plot(V_A, n_A_list, "black")
# plt.plot(V_H, n_H_list, "black")
# plt.plot(V_A_flaps, n_A_list_flaps, "black")
# plt.plot(x_values_AD, y_values_AD, "black", marker="o")
# plt.plot(x_values_FE, y_values_FE, 'black', marker="o")
# plt.plot(x_values_DE, y_values_DE, "black", marker="o")
# plt.plot([point_H[0], point_F[0]], [point_H[1], point_F[1]], "black", marker="o")
# plt.plot([V_max_flaps, V_intersect], [2, 2], "black")
# plt.plot([0, V_D], [1, 1], 'black', linestyle="--")
# plt.plot([V_A_fin, point_A[0]], [0, point_A[1]], 'black', linestyle="--")
# plt.plot([V_S, V_S], [0, 1], "black", linestyle="--")
# plt.plot([point_F[0], V_C], [point_F[1], 0], "black", linestyle="--")
# plt.plot([0,0],[0,0],"black",marker="o")
# plt.text(0-0.015, 0+0.25, "")
# plt.text(point_A[0]-2, point_A[1]+0.2, "A")
# plt.text(point_D[0]-2, point_D[1]+0.25, "D")
# plt.text(point_E[0]+2, point_E[1]+0.25, "E")
# plt.text(point_F[0]-1, point_F[1]-0.25, "F")
# plt.text(point_H[0]-2, point_H[1]-0.25, "H")
# plt.text(V_S-2,0-0.2, r'$V_{S}$')
# plt.text(V_A_fin-2,0-0.2, r'$V_{A}$')
# plt.text(V_C-2,0+0.1, r'$V_{C}$')
# plt.text(V_D-2,0-0.25, r'$V_{D}$')
# plt.xlim(left=0, right=(V_D+20))
# plt.ylim(top=3,bottom=-1.5)
# # plt.xlabel("V")
# ax.set_xlabel("V", weight = "bold")
# ax.xaxis.set_label_coords(0.96, 0.3)
# plt.ylabel("n", weight = "bold")
# ax.spines['bottom'].set_position(('data',0))
# plt.grid()
# plt.show()
#------------------------------------------GUST DIAGRAM DESIGN---------------------------------------------------------
#Constants
CL_alpha = CL_Alpha_Wing
chord = bw / Aw
S = m_mto * g / W_S_design
W_S_oem = m_oem * g / S_w


#Gust Profile
if density == 0:
    U_ref_C = 17.07
    U_ref_D = U_ref_C / 2
elif density == 1:
    U_ref_C = 11.37
    U_ref_D = U_ref_C / 2
else:
    U_ref_C = 12.71
    U_ref_D = U_ref_C / 2


Z_mo = h_cruise
R1 = beta_s_land_fc
R2 = m_zf/m_mto
F_gz = 1 - (Z_mo / 76200)
F_gm = np.sqrt(R2 * np.tan((np.pi * R1)/4))
F_g = 0.5 * (F_gz + F_gm)
H1 = c_mac_w * 12.5
H2 = 107
if H1 < H2:
    H = H2
else:
    H = H1

U_ds = U_ref_C * F_g * (H/107)**(1/6)
U_list = []
H_list = []
for i in np.arange(0,2*H+1,1):
    H_list.append(i)
    U = (U_ds/2)*(1-np.cos((np.pi*i/H)))
    U_list.append(U)

#Design speed for max. gust intensity
mu = (2* W_S_oem) / (rho * c_mac_w * CL_alpha * g)
K_G = (0.88 * mu) / (5.3 + mu)
V_B = V_S * np.sqrt(1 + ((K_G * rho_0 * U_ref_C * V_C * CL_alpha)/(2* W_S_oem)))

V = V_D * np.sqrt(rho_0/rho) # CHANGE BASED ON VELOCITY
lmbda = (2 * W_S_oem) / (CL_alpha * rho * V * g)

U_ref_V = U_ref_D # CHANGE BASED ON VELOCITY
H = 9
dH = 1
max_gust = []
while H <= 107:
    omega = (np.pi * V) / H
    U_ds = U_ref_V * F_g * (H/107)**(1/6)
    t = 0.00005
    dt = 0.005
    ns = []
    nsmax = []
    while t < (2*np.pi)/omega:
        dn = (U_ds / (2*g))*(omega * np.sin(omega * t)+(1/(1+(omega * lmbda)**(-2)))*((np.exp(-t/lmbda)/lmbda)-(np.cos(omega*t)/lmbda)-(omega*np.sin(omega*t))))
        ns.append(dn)
        t = t + dt
        nmax = max(ns)
    nsmax.append(nmax)
    H_list.append(H)
    H = H + dH
    max_gust.append(nmax)

#print("position of max delta n=", np.argmax(max_gust))
# print("delta n =", max_gust[18])



#Velocity values at OEM
V_B_max_dn = 1.86407823
V_C_max_dn = 2.24796708
V_D_max_dn = 1.26166518


#Coordinates for Gust Diagram
pt1 = [0, 1]
pt2 = [V_B, 1+V_B_max_dn]
pt3 = [141.47, 1+V_C_max_dn]
pt4 = [158.8, 1+V_D_max_dn]
pt5 = [158.8, 1-V_D_max_dn]
pt6 = [141.47, 1-V_C_max_dn]
pt7 = [V_B, 1-V_B_max_dn]
pt8 = [0, 1]


x_values = [pt1[0], pt2[0], pt3[0], pt4[0], pt5[0], pt6[0], pt7[0], pt8[0]]
y_values = [pt1[1], pt2[1], pt3[1], pt4[1], pt5[1], pt6[1], pt7[1], pt8[1]]
# print(min(y_values))
# print(max(y_values))
# plt.plot(x_values, y_values, 'bo', linestyle="--")
# plt.axhline(y = 1, color = 'r', linestyle = '-')
# plt.show()

#Find intersections

#------------------------------FINAL COMBINED PLOT---------------------------------------------------------------------
fig, ax = plt.subplots()
plt.rcParams['axes.labelsize'] = 20
plt.axhline(0, color = '0.75')
#plt.tick_params(labelbottom = False)
plt.plot(V_A, n_A_list, "black")
plt.plot(V_H, n_H_list, "black")
plt.plot(V_A_flaps, n_A_list_flaps, "black")
plt.plot([V_max_flaps, V_intersect], [2, 2], "black")
plt.plot(x_values_AD, y_values_AD, "black", marker="o")
plt.plot(x_values_FE, y_values_FE, 'black', marker="o")
plt.plot(x_values_DE, y_values_DE, "black", marker="o")
plt.plot([point_H[0], point_F[0]], [point_H[1], point_F[1]], "black", marker="o")
plt.plot([0, V_S], [1, 1], 'black', linestyle=":")
plt.plot([V_A_fin, point_A[0]], [-1, point_A[1]], 'black', linestyle=":")
plt.plot([V_S, V_S], [-1, 1], "black", linestyle=":")
plt.plot([V_B, V_B], [1+V_B_max_dn, 1-V_B_max_dn], "black", linestyle=":")
plt.plot([point_F[0], V_C], [point_F[1], 1+V_C_max_dn], "black", linestyle=":")
plt.plot([0,0],[0,0],"black",marker="o")
plt.plot(x_values, y_values, 'black', marker="o", linestyle="--")
plt.plot([0, 0], [0, 1], "black", linestyle = '-')
plt.text(0-0.015, 0+0.25, "")
plt.text(V_S-3,0.5, r'$V_{S}$')
plt.text(V_A_fin-3,0.5, r'$V_{A}$')
plt.text(V_C-3,0.5, r'$V_{C}$')
plt.text(V_D-3,0.5, r'$V_{D}$')
plt.text(V_B-3,0.5, r'$V_{B}$')
plt.xlim(left=0, right=(V_D+20))
plt.ylim(top=4,bottom=-1.5)
plt.xlabel("V")
plt.xlabel("Equivalent Airspeed, $V_{EAS}$ m/s", weight = "bold")
ax.xaxis.set_label_coords(0.5, 0)
plt.ylabel("Load Factor, n", weight = "bold")
ax.spines['bottom'].set_position(('data',0))
plt.grid()
ax.set_xlim(xmin=0)
plt.show()
