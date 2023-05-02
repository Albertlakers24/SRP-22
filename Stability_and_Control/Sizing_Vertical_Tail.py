import numpy as np
import math as m
from sympy import symbols, Eq,nsolve
import sympy as sp
from Constants.AircraftGeometry import l_f, d_f_outer, Aw, S_w, bw, taperw, Sweep_quarterchordw, c_mac_w, t_c_ratio_w
from Constants.Aircraft_Geometry_Drawing import center_S_BS,y_T_outboard, zv, y_T, Zw, S_BS, center_S_BS
from Constants.Masses_Locations import m_pldes, xcg_aft_potato, xcg_front_potato,LEMAC, m_mto
from Constants.Aerodynamics import Alpha_DesTakeOff,R_lfus, CL_Alpha_VT, Cl_Alpha_VT_Airfoil, CD0_CR,CL_DesCruise, CL_MaxTakeOff
from Constants.MissionInputs import V_approach, ISA_calculator, landing_critical, dt_land, V_cruise, g, rho_5000
from Constants.FlightPerformance_Propulsion import Power_tot

print("FILE: Vertical_Tail")
"""
VT design checks:
VT.1. Directional stability:
        For sideslip angles: low AR, dorsal fin, large leading edge sweep
        Cn_beta > 0 and Cnr < 0
VT.2. Spin recovery       -> This is a transport aircraft and therefore spin recovery is not necessary

Rudder Design checks:
R.1. Multiple engine inoperative (Directional Trim)
R.2. Cross wind recovery
R.3. Spin recovery        -> This is a transport aircraft and therefore spin recovery is not necessary (??)
"""

CheckVT = 0
CheckR = 0

# Estimated Inputs for VT and Rudder design                                     ITERATIVE PROCESS RESULTS PRESENTED IN THE REPORT
V_v = 0.07                      # VT : volume fraction              [m^3]
l_v = 8.3                       # VT : tail arm                     [m]
Av = 1.5                        # VT : aspect ratio                 [-]
taperv = 0.6                    # VT : taper ratio                  [-]
Sweepv = 26*np.pi/180           # VT : sweep angle                  [rad]
i_v = 0                         # VT : incidence angle              [rad]       -> Symmetric propulsion
Dihedral = 0                    # VT : dihedral                     [rad]       -> Symmetric propulsion
br_bv = 0.9                     # Rudder-to-VT span                 [-]         -> rudder effectiveness, 0.8
Cr_Cv = 0.3                     # Rudder-to-VT chord                [-]         -> graph 12.12 from book
V_stall_TO = np.sqrt(m_mto*g/(S_w*0.5*rho_5000*CL_MaxTakeOff))
V_mincont = 0.8 * V_stall_TO    # Min controllable speed            [m/s]       -> See FAR regulations (estimate 80% of stall speed at take off)
delta_r_max = 30                # Rudder : max allowable deflection [deg]
Cdy = 0.8                       # Aircraft side drag coefficient    [-]         (0.5-0.8)
y_i_rudder = 0                  # Rudder : Inboard location         [m]

# Initial Geometry Calculations for VT
def VT_Geometry():
    Sv = V_v * (S_w * bw) / l_v
    bv = np.sqrt(Av * Sv)                   # VT : Wing span         [m]
    c_mac_v = Sv / bv                       # VT : MAC               [m]
    c_rv = 2 * Sv / ((1 + taperv) * bv)     # VT : Root chord        [m]
    c_tv = taperv * c_rv                    # VT : Tip chord         [m]
    return bv, c_mac_v, c_rv, c_tv, Sv

# Variables with new names
h1 = d_f_outer              # Fuselage : height at h1               [m]
h2 = d_f_outer              # Fuselage : height at h2               [m]
c_mac_v = VT_Geometry()[1]  # VT :  MAC                             [m]
Cvroot = VT_Geometry()[2]   # VT : root chord                       [m]
Cvtip = VT_Geometry()[3]    # VT : tip chord                        [m]
bv = VT_Geometry()[0]       # VT : span                             [m]
Sv = VT_Geometry()[4]       # VT : surface area                     [m^2]

# Intermediate calculations
y_o_rudder = br_bv*bv       # Rudder: Outboard location             [m]
Sidewash_grad = 0.724 +3.06 *((Sv/S_w)/(1+np.cos(Sweep_quarterchordw)))+0.4*(Zw/d_f_outer)+0.009*Aw    # Sidewash gradient   [-]
rho = ISA_calculator(landing_critical, dt_land)[2]                   # Density at critical landing (5000ft)

print('----------------- NEEDED TO CALCULATE CN BETA -----------------')
print("xm/l_f =", xcg_aft_potato/l_f)
print("l_f**2/S_BS =", l_f**2/S_BS)
print("sqrth1/h2 =", np.sqrt(h1/h2))
print("h/wf =", 1)                          # h/wf = 1 in our aircraft
print("Rfus *10^6", R_lfus*10**(-6))
print("bv/2r1=", bv/2*1.75)

print('----------------- NEEDED TO CALCULATE Cnr -----------------')
print("taper ratio = ", taperw)
print("Sweep c/4 = ", Sweep_quarterchordw*180/np.pi, "deg")
print("xbar/c_mac =", 0)        # Overdesigned
print("AR_wing =", Aw)

print("----------------- NEEDED TO CALCULATE Cy_delta_r -----------")
print("Clalpha/Clalpha_theory = ", Cl_Alpha_VT_Airfoil /(2*np.pi))
print("Cf/C = ", Cr_Cv)
print("t/c_VT =", 0.18)
print("eta_i = ", y_i_rudder/bv, "m")
print("eta_o = ", y_o_rudder/bv, "m")
print("taper ratio VT =", taperv)
print("Highest allowable flap deflection = ", delta_r_max, "deg")

print("----Needed to calculate the effective aspect ratio VT ----")
print("taper_v =", taperv)
print("bv/2r1=", bv/2*1.75)
print("x/cr=", 0.25)

# TODO: read graphs for Cn_beta
K_N = 0.00175            # Empirical factor     [-]     Fig 10.28 (p.431)
K_Rl = 1.4               # Factor               [-]     Fig 10.29 (p.432)
kv = 1.0                 # Empirical factor     [-]     Fig 10.12 (p.417)
alpha_TO = Alpha_DesTakeOff*np.pi/180      # Angle of attack at TO  [rad]

# TODO: read graphs for Cn_r
Cnr_CL2 = -0.2           # Fraction             [-]         Fig 10.44 (p. 465)
Cnr_CD0 = -0.3           # Fraction             [rad^-1]    Fig 10.45 (p. 466)

# TODO: read graphs for Cy_delta_r
k_prime =  0.6                                  # Correction factor     [-]         Fig 8.13 (p. 260)
K_b = 0.8-0.1                                   # Flap-span fraction    [-]         Fig 8.51 (p. 292)
Cldelta_over_Cldeltatheory  = 0.8               # Fraction              [-]         Fig 8.15 (p. 262)
Cldeltatheory = 2.8                             # Fraction              [rad-1]     Fig 8.14 (p. 260)

# TODO: read graphs for effective aspect ratio
A_vf_Av = 1.6             # In isolated tail                            [-]         Fig 10.14 (p.420)
A_vhf_Avf = 1.5           # In the presence of fuselage alone           [-]         Fig 10.15 (p.420)
K_vh = 0.9                # Factor for relative size Sh and Sv          [-]         Fig 10.16 (p.422)   -> overdesigned
Av_eff = A_vf_Av*Av*(1+K_vh*A_vhf_Avf-1)

# TODO: Rudder angle of attack effectiveness
tau_r = 0.52              # Rudder angle of attack effectiveness        [-]         (p.682 Mohammed)

def Deriv_Directional_Stability():
    """
    VT.1.       Directional Stability Check
                -> critical at Take Off (?)
                Cn_r    (p.464)
                Cn_beta (p.429)
    """

    Cybeta_v = -kv * CL_Alpha_VT * Sidewash_grad * (Sv / S_w)                               # eq.10.28 (p.418)
    Cn_beta_w = 0
    Cn_beta_f = -57.3*K_N*K_Rl*((S_BS*l_f)/(S_w*bw))
    Cn_beta_v = -Cybeta_v*(l_v*np.cos(alpha_TO)+zv*np.sin(alpha_TO))/bw                           # Determine alpha
    Cn_beta = Cn_beta_w + Cn_beta_f + Cn_beta_v                                             # eq.10.40 (p.429)
    #Cn_beta = 0.09

    Cnr_w = Cnr_CL2 * CL_DesCruise ** 2 + Cnr_CD0 * CD0_CR
    Cnr_v = (2 / bw ** 2) * ((l_v * np.cos(alpha_TO) + zv * np.sin(alpha_TO)) ** 2) * Cybeta_v
    Cn_r = Cnr_v + Cnr_w                                                                    # eq.10.86 (p.464)
    #Cn_r = -0.18
    return Cn_beta, Cn_r,Cybeta_v

Cn_beta = Deriv_Directional_Stability()[0]
Cn_r = Deriv_Directional_Stability()[1]
Cybeta_v = Deriv_Directional_Stability()[2]

# R.1. One Engine Inoperative AND R.2. Cross-Wind Landing calculations
# R.1 most crictical case at TO
def Deriv_Rudder():
    """
    Cydelta_r    (p.493)    [rad^-1]
    Cndelta_r    (p.494)    [rad^-1]
    """
    Cydelta_r = CL_Alpha_VT*(k_prime*K_b)*Cldelta_over_Cldeltatheory*Cldeltatheory*(Sv/S_w)         # eq.10.123 (p.493)
    #Cydelta_r = 0.25

    # todo : fix this (lateral position of propulsion system) - conservative: 5inboard, 8.5outboard
    # Cndelta_r = - Cydelta_r * (l_v * np.cos(alpha_TO) + zv * np.sin(alpha_TO)) / bw                       # eq.10.125 (p.494)
    eta_v = 0.97
    Cndelta_r = -CL_Alpha_VT*V_v*eta_v*(tau_r)*(br_bv)
    #Cndelta_r = -0.09
    return Cydelta_r, Cndelta_r

Cydelta_r = Deriv_Rudder()[0]
Cndelta_r = Deriv_Rudder()[1]

P_engine = Power_tot/4                  # Power per engine (4 motors)       [W]
T_L = P_engine/V_cruise                 # Force by 1 engine                 [N]
Vw =  52.37                             # Maximum cross-wind speed          [m/s]       -> FAR regulations
V_T = np.sqrt(V_approach**2 + Vw**2)    # Total airspeed                    [m/s]
beta = m.atan(Vw/V_approach)            # Side slip angle                   [rad]
Fw = 0.5*rho*(Vw**2)*S_BS*Cdy           # Force generated by cross wind     [N]
N_A = -T_L * y_T_outboard               # Yawing moment for outboard engine inoperative        [N]
Cnzero = 0                              # Cn0                               [-]

print("Side slip angle =", beta*180/np.pi, "deg -> should be > 0 for a mean positive rudder deflection")

"ASYMMETRIC THRUST REQUIREMENT OUTPUT"
delta_r_assym = (T_L*4+T_L*8.5)/(-0.5*rho*(V_mincont**2)*S_w*bw*Cndelta_r)         # Rudder deflection angle        [rad]         ex.12.1 [716 Mohammed]

"CROSS WIND REQUIREMENT OUTPUT"
beta_speed = m.atan(Vw/V_approach)            # [rad]
dc = center_S_BS-xcg_aft_potato               # [m]
Cyzero = 0
Cybeta =-1.5

sigma,delta_r_crosswind = symbols('sigma,delta_r_crosswind')
eq1 = Eq(0.5*rho*(V_T**2)*S_w*bw*(Cnzero+Cn_beta*(beta_speed-sigma) +Cndelta_r*delta_r_crosswind)+Fw*dc*sp.cos(sigma),0)
eq2 = Eq(Fw,0.5*rho*(V_T**2)*S_w*(Cyzero+Cybeta*(beta_speed-sigma)+Cydelta_r*delta_r_crosswind))

sol_dict = nsolve((eq1,eq2),(sigma,delta_r_crosswind), (1,1))
delta_r_crosswind = sol_dict[1]*180/np.pi
crab_angle = sol_dict[0]*180/np.pi
print("crab angle=", crab_angle, "deg")
print("Delta_r_crosswind=", delta_r_crosswind, "deg")

print("----- Directional Stability Requirements ------")
if Cn_beta < 0:
    print("VT does NOT meet directional stability requirement")
    print("Cn_beta = ", Cn_beta, "< 0")
    print("Cn_r = ", Cn_r)
else:
    if Cn_r > 0:
        print("VT does NOT meet directional stability requirement")
        print("Cn_r = ", Cn_r, "> 0")
        print("Cn_beta = ", Cn_beta)
    else:
        print("VT meets directional stability requirement :)")
        print("Cn_beta = ", Cn_beta, "> 0")
        print("Cn_r = ", Cn_r, "< 0")
        CheckVT = 1

print("-----Rudder Requirements-----")
if tau_r > 1:
    print("VT does NOT meet requirement as rudder cannot satisfy directional control/trim req")
    print("As tau_r = ", tau_r , "> 1")
else:
    print("VT meets requirement for rudder directional control/trim requirement: tau_r = ", tau_r, "< 1")
    if abs(delta_r_crosswind) > 30:
        print("Rudder does NOT meet crosswind requirement")
        print("delta_r_crosswind = ", delta_r_crosswind,"> 30 deg")
    else:
        if abs(delta_r_assym) > 30*0.01745:
            print("Rudder does NOT meet asymmetric thrust requirement")
            print("delta_r_assym = ", delta_r_assym*180/np.pi, "> 30 deg")
        else:
            print("Rudder meets crosswind and asymmetric thrust requirements:")
            print("delta_r_crosswind = ", delta_r_crosswind, "deg")
            print("delta_r_asymm = ", delta_r_assym*180/np.pi, "deg")
            print("Rudder effectiveness = ", tau_r)
            CheckR = 1

print("---------CALCULATED DERIVATIVES--------")
print("Cndelta_r =", Cndelta_r)
print("Cydelta_r =", Cydelta_r)
print("Cnbeta =", Cn_beta)
print("Cnr =", Cn_r)
print("A_eff =", Av_eff)
print("Cybeta_v = ",Cybeta_v)

if CheckVT == 1 and CheckR ==1:
    print("---------------------------------------------------------")
    print("Yayy Vertical Tail and Rudder Requirements are met !!")
    print("---- VT Dimensions ----")
    print("Sv = ", Sv, "m^2")
    print("bv = ", bv, "m")
    print("Av =", bv**2/Sv)
    print("Av_eff =", Av_eff)
    print("Cv_r = ", Cvroot, "m")
    print("Cv_t = ", Cvtip, "m")
    print("taperv = ", taperv)
    print("MACV", c_mac_v, "m")

    print("---- Rudder Dimensions ----")
    Crroot = Cr_Cv * Cvroot
    Crtip = Cr_Cv * Cvtip
    br = br_bv * bv
    taper_r =  Crtip/Crroot
    MACr = Crroot * 2/3 * ((1+taper_r + taper_r**2)/ (1 + taper_r))
    Sr = br *  MACr

    print("br =", br)
    print("Croot =", Crroot)            # around 1, Multiple Engine Inoperative (MEI)
    print("Crtip =", Crtip)
    print("taperr", taper_r)
    print("MAC =", MACr)
    print("Sr =", Sr)



"""
# Inputs :  Estimated
lv = l_v
A_v = Av
W_pmax = m_pldes * 9.81
W_to = m_mto * 9.81
rho = 1.225

# Written Variables
Y_e = 4  #m
P_eq = 3138  #Horsepower   
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