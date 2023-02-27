import numpy as np
import math as m
# from Initial_Aircraft_Sizing.Empennage_Design import l_v, Av, Sv
# from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import  m_payload, m_mto
# from Aerodynamic_characteristics.Wing_lift_estimation import Calculate_wingsweep
# from Control_and_Stability.Lateral_Control_derivatives import  Cydelta_r, Cndelta_r


"""
VT design checks:
VT.1. Directional stability:
        For sideslip angles: low AR, dorsal fin, large leading edge sweep
        Cn_beta > 0 and Cnr < 0
VT.2. Spin recovery       -> This is a transport aircraft and therefore spin recovery is not necessary

Rudder Design checks:
R.1. One engine inoperative (Directional Trim)
R.2. Cross wind recovery
R.3. Spin recovery        -> This is a transport aircraft and therefore spin recovery is not necessary (??)
"""


# Imported Variables for VT design
m_payload =1                # Payload mass          [kg]
m_mto =18650                # Maximum TO mass       [kg]
l_f = 23.9                  # Fuselage length       [m]
Sw = 59.9                   # Wing surface area     [m^2]
bw = 26.8                   # Wing span             [m]
D_outer = 3.0               # Outer diameter        [m]
x_cgaft = 12.7              # Cg location aft       [m]
R_lfus = 115629986.920      # Reynolds number fuselage [-]
Sidewash_grad = 1           # Sidewash gradient     [-]
C_L_alpha_v = 1             # CL_alpha VT           [rad^-1]        -> p. 386

# Imported Variables for Cn_beta
h1 = D_outer                                # Height fuselage at h1  [m]
h2 = D_outer                                # Height fuslage at h2   [m]
xm = x_cgaft                                # xcg aft location       [m]
S_BS = 63.89                                # Body side area         [m^2]  TODO: calculate

# Imported Variables for Cn_r
taperw = 1              # Wing taper ratio                  [-]
CLw = 0.63              # Wing lift coefficient             [-]     -> CLdes
Sweep_quarterchord = 1  # Wing sweep at c/4                 [rad]
xac = 1                 # Aerodynamic center location       [m]     -> xac must be aft of cg
c_mac = 2.3             # Wing mean aerodynamic chord       [m]
CD0w =0.021             # Wing zero Drag coefficient of the [-]
ARw = 14                # Wing aspect ratio                 [-]








# Estimated Variable Inputs     -> Based on Initial_Tail_Sizing python file
V_v = 0.08                  # VT volume fraction    [m^3]
l_v = 9.2                   # VT tail arm           [m]
Av = 1.3                    # VT aspect ratio       [-]
Sv = 13.959                 # VT surface area       [m^2]
taperv = 0.5                # VT taper ratio        [-]
i_v = 0                     # VT incidence angle    [rad]       -> Symmetric propulsion
Sweepv = 26*np.pi/180       # VT sweep angle        [rad]
Dihedral = 0                # VT dihedral           [rad]       -> Symmetric propulsion

# Initial Geometry Calculations
def VT_Geometry():
    bv = np.sqrt(Av * Sv)                   # Wing span              [m]
    c_mac_v = Sv / bv                       # MAC                    [m]
    c_rv = 2 * Sv / ((1 + taperv) * bv)     # Root chord             [m]
    c_tv = taperv * c_rv                    # Tip chord              [m]
    return bv, c_mac_v, c_rv, c_tv

# VERTICAL TAIL DESIGN CHECKS
# VT.1. CHECK: Directional Stability
# Cn_beta > 0

print('----------------- NEEDED TO CALCULATE CN BETA -----------------')
print("xm/l_f =", xm/l_f)
print("l_f**2/S_BS =", l_f**2/S_BS)
print("sqrth1/h2 =", np.sqrt(h1/h2))
print("h/wf =", 1)                          # h/wf = 1 in our aircraft
print("Rfus *10^6", R_lfus*10**(-6))

# TODO: read graphs for Cn_beta
K_N = 0.00125            # Empirical factor     [-]     Fig 10.28 (p.431)
K_Rl = 1.945             # Factor               [-]     Fig 10.29 (p.432)
alpha = 4*np.pi/180      # Angle of attack at ? [rad]
kv = 0.95                # Empirical factor     [-]     Fig 10.12 (p.417)

# Intermediate calulations
zv = 3.3026              # Vertical distance body axis to ac VT     [m] TODO: determine
Cybeta_v = -kv*C_L_alpha_v* Sidewash_grad* (Sv/Sw)                  # eq. XX (p. XXX)       TODO: which CLalpha is used?

# Cnr < 0
print('----------------- NEEDED TO CALCULATE Cr -----------------')
print("taper ratio = ", taperw)
print("Sweep c/4 = ", Sweep_quarterchord*180/np.pi, "deg")
print("xbar/c_mac =", (xac -x_cgaft)/c_mac)
print("AR_wing =", ARw)

# TODO: read graphs for C_r
Cnr_CL2 = -0.02          # Fraction             [-]         Fig 10.44 (p. 465)
Cnr_CD0 = -0.3           # Fraction             [rad^-1]    Fig 10.45 (p. 466)

def Deriv_Directional_Stability():
    Cn_beta_w = 0
    Cn_beta_f = -57.3*K_N*K_Rl*((S_BS*l_f)/(Sw*bw))
    Cn_beta_v = -Cybeta_v*(l_v*np.cos(alpha)+zv*np.sin(alpha))/bw
    Cn_beta = Cn_beta_w + Cn_beta_f + Cn_beta_v                                             # XX.XX (p. xxx)
    #Cn_beta = 0.09

    Cnr_w = Cnr_CL2 * CLw ** 2 + Cnr_CD0 * CD0w
    Cnr_v = (2 / bw ** 2) * ((l_v * np.cos(alpha) + zv * np.sin(alpha)) ** 2) * Cybeta_v
    Cn_r = Cnr_v + Cnr_w                                                                    # 10.86 (p.464)
    #Cn_r = -0.18
    return Cn_beta, Cn_r

Cn_beta = Deriv_Directional_Stability()[0]
Cn_r = Deriv_Directional_Stability()[1]

print("----- Directional Stability Requirements ------")
print("Cn_beta =",Cn_beta, "Cn_beta > 0")
print("Cn_r =",Cn_r, "Cn_r < 0")

# RUDDER DESIGN CHECKS
# R.1. One Engine Inoperative
# Inputs : Imported or previously calculated
T_L = 4164                      # Force by main engine          [N]
y_T = 4                         # ?? [m]
CL_alpha_vtail = 2.55           # Wing Lift slope for VT        [rad^-1]
Sv_ini = Sv                     # Initial VT surface area       [m^2]
Cvroot = VT_Geometry()[2]       # VT root chord                 [m]
Cvtip = VT_Geometry()[3]        # VT tip chord                  [m]
bv = VT_Geometry()[0]           # VT span                       [m]
rho = 1.225                     # Density at approach               TODO: to be revised for approach altitude

# Inputs : Estimated
br_bv = 0.7                     # rudder rudder-to-VT span      [-]
V_mincont = 0.8 * 58.8/1.05     # Min controllable speed        [m/s]   -> See FAR regulations (estimate 80% of stall speed)
Cr_Cv = 0.22                    # Rudder-to-VT chord ratio      [-]     -> graph 12.12 from book
eta_v = 0.97                    # VT dynamic pressure ratio     [-]     -> q_v/q_inf    TODO: estimate or calculate?
# delta_r = - np.pi / 6           # Max rudder deflection angle   [rad]

# Intermediate calculations
N_A = -T_L * y_T                # Yawing moment                 [N]

## Needed for the rudder derivative calculations:
# Inputs : Imported or previously calculated
Cl_alpha_vtail = 6.4588         # VT: airfoil lift over alpha slope     [rad^-1]
t_over_c_vtail = 0.15           # VT: Thickness over chord ratio        [-]
y_i_rudder = 0                  # Rudder: Inboard location              [m]
y_o_rudder = br_bv*bv           # Rudder: Outboard location             [m]
delta_r_max = 30                # Rudder: Highest allowable deflection  [deg]

print("----Needed to calculate the Rudder derivatives ----")
print("Clalpha/Clalpha_theory = ", Cl_alpha_vtail/(2*np.pi))
print("Cf/C = ", Cr_Cv)
print("t/c =", t_over_c_vtail)
print("eta_i = ", y_i_rudder/bv, "m")
print("eta_o = ", y_o_rudder/bv, "m")
print("taper ratio VT =", taperv)
print("Highest allowable flap deflection = ", delta_r_max, "deg")

# Inputs : From Graphs
k_prime =  0.675                                # Correction factor     [-]         Fig 8.13 (p. 260)
K_b = 0.7                                       # Flap-span fraction    [-]         Fig 8.51 (p. 292)
Cldelta_over_Cldeltatheory  = 1.0               # Fraction              [-]         Fig 8.15 (p. 262)
Cldeltatheory = 3.5                             # Fraction              [rad-1]     Fig 8.14 (p. 260)

def Deriv_Rudder():             # TODO: to be calculated
    Cydelta_r = C_L_alpha_v*(k_prime*K_b)*Cldelta_over_Cldeltatheory*Cldeltatheory*(Sv/Sw)          # eq. XX (p. XXX)
    Cydelta_r = 0.25

    Cndelta_r = - Cydelta_r * (l_v * np.cos(alpha) + zv * np.sin(alpha)) / bw                       # eq. XX (p. XXX)
    Cndelta_r = -0.09
    return Cydelta_r, Cndelta_r

Cydelta_r = Deriv_Rudder()[0]
Cndelta_r = Deriv_Rudder()[1]

tau_r = Cndelta_r / (br_bv * (( l_v * Sv_ini /(Sw*bw))) * -CL_alpha_vtail * eta_v)      # Rudder angle of attack effectiveness (eq. 12.100)

delta_r_assym = (T_L*y_T)/(-0.5*rho*(V_mincont**2)*Sw*bw*Cndelta_r)                     # Rudder deflection angle        [rad]

# R.2. Cross-Wind Landing
# Inputs : Imported
Vw = 52.37                      # Maximum cross-wind speed      [m/s]   -> FAR regulations
V_approach = 60                 # Approach speed                [m/s]
Ss = S_BS                       # Side area of the aircraft     [m^2]
dc = 1.5                        # Distance center side area to xcg_aft or front (depensd on what's more important)  [m]     TODO: revisit
Cdy = 0.65                      # Aircraft side drag coefficient    [-]     (0.5-0.8)

# Intermediate calculations
V_T = np.sqrt(V_approach**2 + Vw**2)    # Total airspeed        [m/s]
beta = m.atan(Vw/V_approach)            # Side slip angle       [rad]
Fw = 0.5*rho*(Vw**2)*Ss*Cdy             # Force generated by cross wind     [N]
sigma = m.acos(-N_A/(Fw*dc))            # Crab angle            [rad]
Cnzero = 0                              # Cn0                   [-]

delta_r_crosswind = ((N_A/(0.5*rho*V_T**2*Sw*bw)) - Cnzero  -Cn_beta*(beta-sigma))*(1/Cndelta_r)        # Rudder deflection angle       [rad]

if tau_r > 1:
    print("VT needs to be redesigned as rudder cannot satisfy directional control/trim req")
else:
    if abs(delta_r_crosswind) > 30*0.01745:
        print("delta_r due to cross wind too big!")
    else:
        if abs(delta_r_assym) > 30*0.01745:
            print("delta_r due to asymmetric thrust is too big!")
        else:
            print("Rudder/VT sizings meet requirements:")
            print("delta_r_crosswind = ", delta_r_crosswind*180/np.pi, "deg")
            print("delta_r_asymm = ", delta_r_assym*180/np.pi, "deg")
            print("Rudder effectiveness = ", tau_r)


print("---- VT Dimensions ----")
# print("taperv", taper_v)
# MACv = Cvroot * 2/3 * ((1+taper_v + taper_v**2)/ (1 + taper_v))
# print("MACV", MACv)
# print(taper_v)

print("---- Rudder Dimensions ----")
Crroot = Cr_Cv * Cvroot
Crtip = Cr_Cv * Cvtip
br = br_bv * bv
taper_r =  Crtip/Crroot
MACr = Crroot * 2/3 * ((1+taper_r + taper_r**2)/ (1 + taper_r))
Sr = br *  MACr
#Sweep_hinge = Calculate_wingsweep(26,0, 3.783, 0.5)    TODO: to be fixed

print("br =", br)
print("Croot =",Crroot)
print("Crtip =",  Crtip)
print("taperr", taper_r)
print("MAC =", MACr)
print("Sr =", Sr)
# print("Sw h", Sweep_hinge)


"""
# Inputs :  Estimated
lv = l_v
A_v = Av
W_pmax = m_payload * 9.81
W_to = m_mto * 9.81
rho = 1.225

# Written Variables
Y_e = 4  #m
P_eq = 3138  #Horsepower   #To be imported from Tims code
C_lmax = 2.1
k_dr = 1.1  #Obtained from graph
k_v = 1.1   #Preset value for our configuration

print(A_v)
C_Y_dr = Cydelta_r
#Calculation of Lateral force needed by the tail
Y_v = - N_A / l_v
print("lv", l_v)
#Calculation of Sv from Y_v

Sv_new = 2 * Y_v / (rho * C_Y_dr * -delta_r  * V**2)

print("S", Sv_new)
print("CLEAR NOW________________________________________")

print("Svini", Sv_ini)
print("bvini", bv)
print("MACv", MACv)
print("Cvr", Cvroot)
print("Cvt", Cvtip)
print("Sr", Sr)
print("br", br)
print("MACr", MACr)
print("Crr", Crroot)
print("Crt", Crtip)

print("NEW TAIL VALUES")
print("Svfin", Sv_new)

AR_v = bv**2 / Sv_ini
bvnew = (AR_v * Sv_new)**0.5
MACv_new = bvnew / AR_v

MACv_ratio = MACv_new / MACv

Cvroot_new = MACv_ratio * Cvroot
Cvtip_new = MACv_ratio * Cvtip

print("bvnew", bvnew)
print("MACvnew", MACv_new)
print("Cvrnew", Cvroot_new)
print("Cvtnew", Cvtip_new)

sweep = Calculate_wingsweep(26,0,3.92,0.5)
print("SWEEP", sweep)

print(l_v)
"""