import numpy as np
import math as m
from Constants.AircraftGeometry import l_f, d_f_outer, Aw, S_w, bw, taperw, Sweep_quarterchordw, c_mac_w, t_c_ratio_w
from Constants.Aircraft_Geometry_Drawing import zv, y_T, Zw, S_BS, center_S_BS
from Constants.Masses_Locations import m_mto, m_pldes
from Constants.Aerodynamics import R_lfus, CL_Alpha_VT, Cl_Alpha_VT_Airfoil
from Constants.MissionInputs import V_approach, ISA_calculator, landing_critical, dt_land

# Imported Variables :  Aerodynamics
CLw = 0.63                  # Wing : lift coefficient               [-]         -> CLdes        todo: cruise? TO? Land?
CD0w =0.021                 # Wing : zero Drag coefficient          [-]                         todo: cruise? deflection flap?
# Imported Variables : Performance
T_L = 4164                  # Force by 1 engine                     [N]
eta_v = 0.97                # VT dynamic pressure ratio             [-]         -> q_v/q_inf    TODO: estimate or calculate?
Vw = 52.37                  # Maximum cross-wind speed              [m/s]   -> FAR regulations
#Imported Variables: Masses & Locations
m_payload = m_pldes         # Payload mass   ?                      [kg]    todo: ask if correct to Albert
x_cgaft = 12.7              # Aircraft aft cg                       [m]     todo: xcg potato or normal??
x_cgfront = 11.7            # Airfract front cg                     [m]
# Imported Variables : Masses and Locations
xac = 13                    # AC location                           [m]         -> xac must be aft of cg    TODO: to be calculated


print("FILE: Vertical_Tail")
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

CheckVT = 0
CheckR = 0

# Estimated Inputs for VT and Rudder design                                     ITERATIVE PROCESS RESULTS PRESENTED IN THE REPORT
V_v = 0.08                      # VT : volume fraction              [m^3]
l_v = 9.2                       # VT : tail arm                     [m]
Av = 1.3                        # VT : aspect ratio                 [-]
taperv = 0.5                    # VT : taper ratio                  [-]
Sweepv = 26*np.pi/180           # VT : sweep angle                  [rad]
i_v = 0                         # VT : incidence angle              [rad]       -> Symmetric propulsion
Dihedral = 0                    # VT : dihedral                     [rad]       -> Symmetric propulsion
br_bv = 0.7                     # Rudder-to-VT span                 [-]
Cr_Cv = 0.22                    # Rudder-to-VT chord                [-]         -> graph 12.12 from book
V_mincont = 0.8 * 58.8/1.05     # Min controllable speed            [m/s]       -> See FAR regulations (estimate 80% of stall speed)    todo: or at TO speed??
delta_r_max = 30                # Rudder : max allowable deflection [deg]
Cdy = 0.65                      # Aircraft side drag coefficient    [-]         (0.5-0.8)
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
rho = ISA_calculator(landing_critical, dt_land)                     # Density at critical landing (5000ft)

print('----------------- NEEDED TO CALCULATE CN BETA -----------------')
print("xm/l_f =", x_cgaft/l_f)
print("l_f**2/S_BS =", l_f**2/S_BS)
print("sqrth1/h2 =", np.sqrt(h1/h2))
print("h/wf =", 1)                          # h/wf = 1 in our aircraft
print("Rfus *10^6", R_lfus*10**(-6))

print('----------------- NEEDED TO CALCULATE Cr -----------------')
print("taper ratio = ", taperw)
print("Sweep c/4 = ", Sweep_quarterchordw*180/np.pi, "deg")
print("xbar/c_mac =", (xac -x_cgaft)/c_mac_w)
print("AR_wing =", Aw)

print("----Needed to calculate the Rudder derivatives ----")
print("Clalpha/Clalpha_theory = ", Cl_alpha_vtail/(2*np.pi))
print("Cf/C = ", Cr_Cv)
print("t/c =", t_c_ratio_w)
print("eta_i = ", y_i_rudder/bv, "m")
print("eta_o = ", y_o_rudder/bv, "m")
print("taper ratio VT =", taperv)
print("Highest allowable flap deflection = ", delta_r_max, "deg")

# TODO: read graphs for Cn_beta
K_N = 0.00125            # Empirical factor     [-]     Fig 10.28 (p.431)
K_Rl = 1.945             # Factor               [-]     Fig 10.29 (p.432)
alpha = 4*np.pi/180      # Angle of attack at ? [rad]
kv = 0.95                # Empirical factor     [-]     Fig 10.12 (p.417)

# TODO: read graphs for Cn_r
Cnr_CL2 = -0.02          # Fraction             [-]         Fig 10.44 (p. 465)
Cnr_CD0 = -0.3           # Fraction             [rad^-1]    Fig 10.45 (p. 466)

# TODO: read graphs for Cy_delta_r
k_prime =  0.675                                # Correction factor     [-]         Fig 8.13 (p. 260)
K_b = 0.7                                       # Flap-span fraction    [-]         Fig 8.51 (p. 292)
Cldelta_over_Cldeltatheory  = 1.0               # Fraction              [-]         Fig 8.15 (p. 262)
Cldeltatheory = 3.5                             # Fraction              [rad-1]     Fig 8.14 (p. 260)

def Deriv_Directional_Stability():
    """
    VT.1.       Directional Stability Check
    :return: Cn_beta and Cn_r
    """

    Cybeta_v = -kv * CL_alpha_vtail * Sidewash_grad * (Sv / S_w)                                # eq. XX (p. XXX)
    Cn_beta_w = 0
    Cn_beta_f = -57.3*K_N*K_Rl*((S_BS*l_f)/(S_w*bw))
    Cn_beta_v = -Cybeta_v*(l_v*np.cos(alpha)+zv*np.sin(alpha))/bw                           # Determine alpha
    Cn_beta = Cn_beta_w + Cn_beta_f + Cn_beta_v                                             # XX.XX (p. xxx)
    #Cn_beta = 0.09

    Cnr_w = Cnr_CL2 * CLw ** 2 + Cnr_CD0 * CD0w
    Cnr_v = (2 / bw ** 2) * ((l_v * np.cos(alpha) + zv * np.sin(alpha)) ** 2) * Cybeta_v
    Cn_r = Cnr_v + Cnr_w                                                                    # 10.86 (p.464)
    #Cn_r = -0.18
    return Cn_beta, Cn_r

Cn_beta = Deriv_Directional_Stability()[0]
Cn_r = Deriv_Directional_Stability()[1]

# R.1. One Engine Inoperative AND R.2. Cross-Wind Landing calculations
def Deriv_Rudder():
    Cydelta_r = CL_alpha_vtail*(k_prime*K_b)*Cldelta_over_Cldeltatheory*Cldeltatheory*(Sv/S_w)       # eq. XX (p. XXX)       TODO: revisit
    #Cydelta_r = 0.25

    Cndelta_r = - Cydelta_r * (l_v * np.cos(alpha) + zv * np.sin(alpha)) / bw                       # eq. XX (p. XXX)       TODO: check for what alpha?
    #Cndelta_r = -0.09
    return Cydelta_r, Cndelta_r

Cydelta_r = Deriv_Rudder()[0]
Cndelta_r = Deriv_Rudder()[1]
print("Cndeltar =", Cndelta_r)
print("Cydelta_r =", Cydelta_r)


tau_r = Cndelta_r / (br_bv * (( l_v * Sv /(S_w*bw))) * -CL_alpha_vtail * eta_v)      # Rudder angle of attack effectiveness (eq. 12.100)     [-]

V_T = np.sqrt(V_approach**2 + Vw**2)    # Total airspeed                    [m/s]
beta = m.atan(Vw/V_approach)            # Side slip angle                   [rad]
Fw = 0.5*rho*(Vw**2)*S_BS*Cdy           # Force generated by cross wind     [N]
N_A = -T_L * y_T                        # Yawing moment for outboard engine inoperative        [N]
Cnzero = 0                              # Cn0                               [-]

print("Side slip angle =", beta)            #TODO: positive sideslip angle should mean positive rudder deflection!

def crab_angle(xcg):
    """crab angle calculations [rad] ;
    :param dc :  distance center side area to xcg aft or front depending on what is critical [m]"""
    dc = xcg - center_S_BS
    sigma = m.acos(-N_A / (Fw * dc))
    return sigma

"ASYMMETRIC THRUST REQUIREMENT OUTPUT"
delta_r_assym = (-T_L*y_T)/(-0.5*rho*(V_mincont**2)*S_w*bw*Cndelta_r)         # Rudder deflection angle        [rad] todo: check sign, as otherwise - can be left out.

"CROSS WIND REQUIREMENT OUTPUT"
delta_r_crosswind_front = ((N_A/(0.5*rho*V_T**2*S_w*bw)) - Cnzero  -Cn_beta*(beta-crab_angle(x_cgfront)))*(1/Cndelta_r)    # Rudder deflection angle for front cg      [rad]
delta_r_crosswind_aft = ((N_A/(0.5*rho*V_T**2*S_w*bw)) - Cnzero  -Cn_beta*(beta-crab_angle(x_cgaft)))*(1/Cndelta_r)        # Rudder deflection angle for aft cg        [rad]
if abs(delta_r_crosswind_aft) > abs(delta_r_crosswind_front):
    delta_r_crosswind = delta_r_crosswind_aft
else:
    delta_r_crosswind = delta_r_crosswind_front

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
    print("tau_r = ", tau_r)
else:
    print("VT meets requirement for rudder directional control/trim requirement: tau_r = ", tau_r, "< 1")
    if abs(delta_r_crosswind) > 30*0.01745:
        print("Rudder does NOT meet crosswind requirement")
        print("delta_r_crosswind = ", delta_r_crosswind,"> 30 deg")
    else:
        if abs(delta_r_assym) > 30*0.01745:
            print("Rudder does NOT meet asymmetric thrust requirement")
            print("delta_r_assym = ", delta_r_assym, "> 30 deg")
        else:
            print("Rudder meets crosswind and asymmetric thrust requirements:")
            print("delta_r_crosswind = ", delta_r_crosswind*180/np.pi, "deg")
            print("delta_r_asymm = ", delta_r_assym*180/np.pi, "deg")
            print("Rudder effectiveness = ", tau_r)
            CheckR = 1

if CheckVT == 1 and CheckR ==1:
    print("---------------------------------------------------------")
    print("Yayy Vertical Tail and Rudder Requirements are met !!")
    print("---- VT Dimensions ----")
    print("Sv = ", Sv, "m^2")
    print("bv = ", bv, "m")
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
    print("Croot =", Crroot)
    print("Crtip =", Crtip)
    print("taperr", taper_r)
    print("MAC =", MACr)
    print("Sr =", Sr)



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


# Imported Variables : Aircraft Geometry
# l_f = 23.9                  # Fuselage : length                     [m]
# d_f_outer = 3.0             # Fuselage : outer diameter             [m]
# Sw = 59.9                   # Wing : surface area                   [m^2]
# bw = 26.8                   # Wing : span                           [m]
# taperw = 1                  # Wing : taper ratio                    [-]
# Sweep_quarterchordw = 0     # Wing : sweep at c/4                   [rad]
# c_mac_w = 2.3               # Wing : MAC                            [m]
# Aw = 14                     # Wing : aspect ratio                   [-]
# t_over_c_vtail = 0.15       # VT : Thickness over chord ratio       [-]
# Imported Variables : Aircraft Geometry Drawing
# zv = 3.3026                 # Vertical distance body axis to ac VT  [m]
# y_T = 4                     # Distance center line fuselage-engine  [m]
# Zw = -1.84                  # Distance center line fuselage-wing    [m]         -> Negative for high wing   (p. 416)
# S_BS = 63.89                # Body side area                        [m^2]
# center_S_BS = 11            # Center location Side area aircraft    [m]
