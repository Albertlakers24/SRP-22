import numpy as np
import math as m
import matplotlib.pyplot as plt
from Constants.Empennage_LandingGear import lh, bh, Sh, lambdahalf_h, c_rh, c_th, A_h, Vh_V, t_over_c_HT, Sweep_quarter_chord_HT
from Constants.AircraftGeometry import S_w, c_mac_w, Aw, Sweep_quarterchordw
from Constants.Masses_Locations import m_mto, LEMAC, xcg_aft_potato
from Constants.MissionInputs import g, ISA_calculator, h_cruise, dt_cruise, M_cruise, V_cruise, V_approach
from Constants.Aerodynamics import CL_Alpha_Wing, downwash, Cl_Alpha_HT_Airfoil, CL_Alpha_HT, Cmac, R_HT
from Constants.Stability_Control import Cmalpha

print("FILE: Horizontal Tail Sizing")


# To be found still :( todo still to be found
AlphaCL0_CR = 0.16                                      # Angle of attack at CL0    [rad]

# Variables for this file:
eta = 0.95                                              #                           [-]
CNh_delta = 0.024*180/np.pi                             # Normal force gradient     [rad^-1]    Graph (pg. 317 FD reader)
i_h = -2*np.pi/180                                      # Incidence angle           [rad]       To minimize parasite drag todo how determined?
x_ac_w_MAC = 0.25                                       # Position ac of the wing   [MAC]

# Intermediate Equations
MTOW = m_mto*g                                          # Maximum to weight         [N]
rho = ISA_calculator(h = h_cruise, dt=dt_cruise)[2]     # Density at cruise         [kg/m^3]
C_N_alpha_h = CL_Alpha_HT                               # Lift curve gradient HT    [rad^-1]

# Graph arange velocities and angle of attack
V_tailload = np.arange(V_cruise-V_cruise*0.7,V_cruise+V_cruise*0.4,0.01)         # Velocity range for tailload      [m/s]
V_controlforce = np.arange(0,250,0.1)                                            # Velocity range for control force [m/s]
V_trim = np.arange(V_cruise-V_cruise*0.7,V_cruise+V_cruise*0.4,0.01)             # Velocity range for trim          [m/s]
alpha_trim = np.arange(-10,10, 0.1)                                              # Alpha range                      [deg]

# Design Choices  TODO: to be chosen :)
deriv_deltae_se = 2.2                                   #                           [rad/m]         normal values: 2.2-26 rad/m (1.25-1.5deg/cm) --- pg 400 FD reader
delta_te = -11*np.pi/180                                # For stable stick free aircraft: delta_te > delta_te0          [rad]            Initial: -11 deg
cf_over_c = 0.3                                         # Flap chord over chord for elevator    [-]
ct_over_c = 0.2                                         # Tab chord over chord for elevator     [-]
outboard = bh/2                                         # Outboard location         [m]
inboard = 0.3                                           # Inboard location          [m]
cb_over_cf = 0.3                                        # Hinge line position       [m]
outboard_tab = 2.5                                      # Outboard location tab     [m]
inboard_tab = 0.5                                       # Inboard location tab      [m]

# Calculations
Vtrim = V_approach                                      # Trim speed                [m/s]
cb_over_c = cb_over_cf*cf_over_c                        # Wing parameter, cb/c      [-]
cbprime_cprime = cb_over_cf*cf_over_c                   # Wing parameter, cb'/c'    [-]
cfprime_cprime = cb_over_c+cf_over_c                    # Wing paramter, cf'/c'     [-]

# Supporting equations sizing elevator
angle = m.atan((c_rh-c_th)/(bh/2))                      # Angle                                 [rad]
root_H = c_rh-inboard*np.tan(angle)                     # Chord of H at start elevator          [m]
Cr_e = root_H*cf_over_c                                 # Root chord elevator                   [m]
MACe = Cr_e*0.75                                        # MAC elevator                          [m]             Estimate??? TODO: revisit
y = np.tan(angle)*(bh/2 - outboard)                     # Difference change                     [m]
Ct_e = (c_th+y)*cf_over_c                               # Tip chord elevator                    [m]
Se = Ct_e*(outboard-inboard)+0.5*(Cr_e-Ct_e)*(outboard-inboard)     # Elevator surface area     [m^2]

print("bh=", bh)
print("Elevator Sizing:")
print("Inboard location=",      inboard)
print("Outboard location=",     outboard)
print("Cr_e=",                  Cr_e)
print("Ct_e=",                  Ct_e)
print("Se=",                    Se, "1 elevator")
print("MACe=", MACe)

# Supporting equations sizing trim tab
root_H_tab = c_rh-inboard_tab*np.tan(angle)             # Chord of H at start trimtab           [m]
Cr_tab = root_H_tab*ct_over_c                           # Root chord trim tab                   [m]
y = np.tan(angle)*(bh/2 - outboard_tab)                 # Difference change                     [m]
Ct_tab = (c_th+y)*ct_over_c                             # Tip chord elevator                    [m]
S_tab = Ct_tab*(outboard_tab-inboard_tab)+0.5*(Cr_tab-Ct_tab)*(outboard_tab-inboard_tab)        # Elevator surface area     [m^2]

print("Trim Tab Sizing:")
print("Inboard location tab=", inboard_tab)
print("outboard location tab=", outboard_tab)
print("Cr_tab=", Cr_tab)
print("Ct_tab=", Ct_tab)
print("S_tab=", S_tab, "1 trim tab")

# Supporting Values from the Graphs for elevator
cla_over_cla_theory_HT = 0.76                           # Graph 10.64a (p. 501)
c_accent_alpha_over_chalpha_theory_HT = 0.3             # Graph 10.63a (p. 500)
chalpha_theory = -0.48                                  # Graph 10.63b (p. 500)
Ch_alpha_bal_over_Chalpha = 0.85                        # Graph 10.65a (p.503)
c_accent_delta_over_chdelta_theory_HT = 0.75            # Graph 10.69a (p. 507)
chdelta_theory_HT = -0.88                               # Graph 10.69b (p. 507)
Ch_delta_bal_over_Chdelta = 0.75                        # Graph 10.71 (p. 509)
tc_over_two_cf = 0.2485                                 # tc = thickness chord at hingeline and cf = chord of flap (p. 509)

print("-------Needed for Graphs for Chdelta and Chalpha-------")
print("cf/c=", cf_over_c)
print("ct/c=", ct_over_c)
print("ct/cf=", ct_over_c*(1/cf_over_c))
print("Needed for Graphs for Chdelta and Chalpha")
print("Reynolds Number for HT=", R_HT)
print("tan(0.5*phi_te)=", t_over_c_HT, "=t/c")
print("cla_over_cla_theory_HT=",cla_over_cla_theory_HT)
print("Determine tc/2cf -> for the balance ratio")
print("Balance Ratio=", np.sqrt((cb_over_cf)**2-(tc_over_two_cf)**2))

def HingemomentAF_Coefficients():
    """ ELEVATOR
    :return:
    Chdelta (rad^-1)    p.506
    Chalpha (rad^-1)    p.499
    """
    c_accent_h_alpha= c_accent_alpha_over_chalpha_theory_HT*chalpha_theory
    Ch_alpha_bal = c_accent_h_alpha*Ch_alpha_bal_over_Chalpha
    ChAF_alpha = Ch_alpha_bal/((1-M_cruise**2)**0.5)

    c_accent_h_delta = c_accent_delta_over_chdelta_theory_HT*chdelta_theory_HT
    ch_delta_bal = c_accent_h_delta*Ch_delta_bal_over_Chdelta
    ChAF_delta = ch_delta_bal/((1-M_cruise**2)**0.5)
    #Ch_alpha = -0.02
    RatioAF = ChAF_alpha/ChAF_delta
    return ChAF_delta, ChAF_alpha, RatioAF

print("-------Needed for Graphs for Chdelta_t-------")
print("ct/cf_trimtab =", "idk")
print("cf/c =","idk")

Chdeltat_cl_delta = -0.013*180*(1/np.pi)                  # Graph 10.72  (p.511)        [rad^-1]
Chcl_deltat_delta = -0.11                                 #
cla_over_cla_theory_trim = 0.76                           # Graph 10.64a (p.501)
c_accent_alpha_over_chalpha_theory_trim = 0.21            # Graph 10.63a (p.500)
chalpha_theory_trim = -0.38                               # Graph 10.63b (p.500)
Ch_alpha_bal_over_Chalpha_trim = 0.85                     # Graph 10.65a (p.503)
alpha_delta_cl_delta = -0.5                               # Change in angle of attack due to change in tab deflection                         [-]

def HingemomentAF_Coefficients():
    """ TRIMTAB
    :return:
    Chdelta_t (rad^-1)      p.510
    Chalpha_t (rad^-1)      p.
    """
    c_accent_h_alpha= c_accent_alpha_over_chalpha_theory_HT*chalpha_theory
    Ch_alpha_bal = c_accent_h_alpha*Ch_alpha_bal_over_Chalpha
    ChAF_alpha = Ch_alpha_bal/((1-M_cruise**2)**0.5)

    Chdelta_t_AF = Chdeltat_cl_delta-Chcl_deltat_delta*Cl_Alpha_HT_Airfoil*alpha_delta_cl_delta
    RatioAF = ChAF_alpha/Chdelta_t_AF
    return Chdelta_t_AF, ChAF_alpha, RatioAF

Sweep_quarter_h = Sweep_quarter_chord_HT                       # Sweep at quarter chord Htail  [rad]

# From Graphs
DeltaCha_over_clalphaBK = 0.01              # Graph 10.77a (p. 515) at A=4
B2 = 1.25                                   # Graph 10.77c (p. 515)
Eta_i = inboard/(bh/2)                      # Half span inboard location [-]
Eta_o = outboard/(bh/2)                     # Half span outboard location [-]
Kalpha_i = 1.1
Kalpha_o = 4.2
ch_alpha_M = HingemomentAF_Coefficients()[1]
ch_delta_M = HingemomentAF_Coefficients()[0]
DeltaChd_cldBKd = 0.016                     # Graph 10.78a (p. 517) at A=4
cl_delta = 5.4                              # Graph 8.14 (p. 260)
Kdelta_o = 4.2                              # Graph 10.78c (p. 517)             normal 3.5
Kdelta_i = 1.05                             # Graph 10.78c (p. 517)
Sweep_hl = 0.2                              # Sweep at the hingline
alpha_d = 0.51                              # Graph 8.17 (p. 262) at delta_f = 35deg

print("Needed for 3D hingemoment calculations")
print("cb'/cf'=", cbprime_cprime)
print("cf'/c'=", cfprime_cprime)
print("cf/c=", cf_over_c)
print("Eta Inboard Location", Eta_i)
print("Eta Outboard Location,", Eta_o)

def Hingemoment_Coefficients():
    Kalpha = Kalpha_i*(1-Eta_i)-Kalpha_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaC_h_alpha = DeltaCha_over_clalphaBK*(Cl_Alpha_HT_Airfoil*B2*Kalpha*np.cos(Sweep_quarter_h))
    Ch_alpha = ((A_h*np.cos(Sweep_quarter_h))/(A_h+2*np.cos(Sweep_quarter_h)))*ch_alpha_M + DeltaC_h_alpha

    Kdelta = Kdelta_i*(1-Eta_i)-Kdelta_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaCh_delta =DeltaChd_cldBKd*(cl_delta*B2*Kdelta*np.cos(Sweep_quarter_h)*np.cos(Sweep_hl))
    Ch_delta = np.cos(Sweep_quarter_h)*np.cos(Sweep_hl)*(ch_delta_M +alpha_d*ch_alpha_M)*((2*np.cos(Sweep_quarter_h))/(A_h+2*np.cos(Sweep_quarter_h)))+DeltaCh_delta

    Ch_delta_horn = Ch_delta*0.26                                    # Horn Effect
    Ratio = Ch_alpha / Ch_delta_horn
    return Ch_delta_horn, Ch_alpha, Ratio

Ch_delta = Hingemoment_Coefficients()[0]
Ch_alpha = Hingemoment_Coefficients()[1]
print("Hinge line coefficients",Ch_delta, Ch_alpha, Hingemoment_Coefficients()[2])
Chalpha_Chdelta = Hingemoment_Coefficients()[2]
if Hingemoment_Coefficients()[0]<0:
    print("Chdelta<0 statisfied", "Chdelta=", Ch_delta)
else:
    print("Chdelta>0 SHOULD BE CHANGED!!")

Eta_i_trim = inboard_tab/(bh/2)
Eta_o_trim = outboard_tab/(bh/2)
print("Eta_i_trim=", Eta_i_trim)
print("Eta_o_trim=", Eta_o_trim)
ch_alpha_M_trim = HingemomentAF_Coefficients()[1]*0.5       # ch_alpha_M for airfoil for trimtab
ch_delta_M_trim = HingemomentAF_Coefficients()[0]*0.2
A_trim = ((outboard_tab-inboard_tab)**2)/S_tab
print("Atrim=", A_trim)
DeltaChd_cldBKd_trim = 0.016                     # Graph 10.78a (p. 517) at A=4
cl_delta_trim = 2.5                              # Graph 8.14 (p. 260)
Kdelta_o_trim = 2.5                              # Graph 10.78c (p. 517)
Kdelta_i_trim = 1.25                             # Graph 10.78c (p. 517)
Sweep_hl_trim = 0.2                              # Sweep at the hingline of trim tab                         [rad]
alpha_d_trim = delta_te                          # Graph 8.17 (p. 262) at delta_f = 35deg of trim tab        [rad]
Sweep_quarter_h_trim = 0.1                       # Quarter chord sweep for trim tab

def Chdelta_t():
    Kdelta = Kdelta_i_trim * (1 - Eta_i_trim) - Kdelta_o_trim * ((1 - Eta_o_trim) / (Eta_o_trim - Eta_i_trim))
    DeltaCh_delta = DeltaChd_cldBKd_trim * (cl_delta_trim * B2 * Kdelta * np.cos(Sweep_quarter_h_trim) * np.cos(Sweep_hl_trim))
    Chdelta_t = np.cos(Sweep_quarter_h_trim) * np.cos(Sweep_hl_trim) * (ch_delta_M_trim + alpha_d_trim * ch_alpha_M_trim) * ((2 * np.cos(Sweep_quarter_h_trim)) / (A_trim + 2 * np.cos(Sweep_quarter_h_trim))) + DeltaCh_delta
    return Chdelta_t

print("Chdelta_t=", Chdelta_t())

# Supporting Equations
CNhalpha_free= C_N_alpha_h - CNh_delta*Chalpha_Chdelta                                                          # Normal force gradient eq. 7.5 Sam1                [rad^-1]
x_ac_w_meters = c_mac_w*x_ac_w_MAC+LEMAC                                                                        # Location aerodynamic center wing                  [m]
xnfree = ((CNhalpha_free/CL_Alpha_Wing)*(1-downwash)*(Vh_V**2)*((Sh*lh)/(S_w*c_mac_w))*c_mac_w) + x_ac_w_meters # Location neutral point stick free eq. 7.7 Sam1    [m]         normal value 0.483*MAC
Cmdelta = -CNh_delta*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))                                                            # eq. 5.21 Sam1                                     [rad^-1]

if xcg_aft_potato <xnfree:
    print("Aircraft is statically stable as xcg<xnfree")
    print("xcg_aft=",xcg_aft_potato )
    print("xnfree=", xnfree)
else:
    print("Aircraft is NOT statically stable as xcg>xnfree")
    print("xcg_aft=", xcg_aft_potato )
    print("xnfree=", xnfree)

def Cm0():
    """ CHECKED
    Cm0 is the Cm at C_N = 0
    :param Cmac: (-)
    :param CN_alpha: for the horizontal tail plane (rad^-1)
    :param AlphaCL0_CR: when CL=0 (rad)
    :param i_h: incidence angle to minimize the parasite drag (rad)
    :return: Cm0 (-)
    """
    Cm0 = Cmac - C_N_alpha_h*(AlphaCL0_CR+i_h)*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))
    return Cm0

print("Where should the suspension point be for stability?")
print("Cmac=",Cmac)

print("Cm0=", Cm0())
if Cm0() > 0:
    print("Cm0 is positive!! Should be negative!")

def Cmdelta_e():
    """ CHECKED
    Control derivative due to elevator deflection equations from FD reader p. 196
    :param CNh_delta: (rad^-1)
    :return: Cmdelta_e: (rad^-1)    normal value around -1.0 to -1.5
    """
    Cmdelta_e = -CNh_delta*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))
    return Cmdelta_e

def delta_e(V):
    """ CHECKED
    Elevator deflection per velocity equation from slide 22 lecture 4 (AE3212)
    :param Cm_delta_e:(rad^-1)
    :param Cm0: (-)
    :param Cm_alpha: (rad^-1)
    :param C_N_alpha (rad^-1)
    :param V: velocity (m/s)
    :return: delta_e: elevator deflection (deg)
    """
    delta_e = -(1/Cmdelta_e())*(Cm0() + (Cmalpha/CL_Alpha_Wing)*(MTOW/(0.5*S_w*rho*V**2)))
    delta_e_degree = delta_e*180/np.pi
    return delta_e_degree

def delta_e_alpha(alpha):
    """ CHECKED
    Elevator deflection per angle of attack equation from slide 24 lecture 4 (AE3212)
    :param Cmdelta_e: (rad^-1)
    :param Cm0: (-)
    :param Cmalpha: (rad^-1)
    :param alpha and AlphaCL0_CR: angle of attack (rad)
    :return: delta_e: elevator deflection (deg)
    """
    alpha_rad = alpha*np.pi/180
    delta_e = -(1/Cmdelta_e())*(Cm0() + (Cmalpha*(alpha_rad-AlphaCL0_CR)))
    delta_e_degree = delta_e*180/np.pi
    return delta_e_degree

def TailLoad(V):
    """ CHECKED
    eq. Lecture 3 slide 33
    Tail Load Diagram, used for rotational and vertical stability,
    important for the structural department
    :param Cmac: (-)
    :param MTOW: maximum take off weight (N)
    :param xcg and xw: locations (m)
    :return: Tail Load (N)
    """
    Nh = (1/lh)*(Cmac*0.5*rho*V**2*S_w*c_mac_w + MTOW*(xcg_aft_potato -x_ac_w_MAC))
    return Nh

def CNh():
    """ CHECKED
    normal force equation
    :return: CNh (-)
    """
    CNh = TailLoad(V = V_cruise)/(0.5*rho*V_cruise**2*S_w)
    return CNh

def trimtab_0():
    """ CHECKED
    :return: delta_te0  (rad)
    """
    delta_te0 = ((Ch_delta /Cmdelta)*Cmac) + ((Ch_delta/CNh_delta)*CNhalpha_free*(AlphaCL0_CR+i_h)) / Chdelta_t()
    return delta_te0

print("delta_t0=", trimtab_0(), "rad")
print("delta_te=", delta_te, "rad")

def ControlForce(V, delta_te):
    """ CHECKED
    Control Curve, used for elevator control force stability (dFe/dV)Fe=0 >0
    :param: Se; Use surface area of 1 elevator (m^2) - same when you have St
    :return: Tail Load (N)
    """
    F_velocity_independent = (MTOW/S_w)*(Ch_delta/Cmdelta_e())*((xcg_aft_potato -xnfree)/c_mac_w)
    F_velocity_dependent= 0.5*rho*V**2*Chdelta_t()*(delta_te-trimtab_0())
    a = deriv_deltae_se*Se*MACe*(Vh_V)**2
    Fe = a*(F_velocity_independent - F_velocity_dependent)
    return Fe

def deriv_controlForce():
    deriv_trim = -2*deriv_deltae_se*Se*MACe*(Vh_V**2)*(MTOW/S_w)*(Ch_delta/Cmdelta_e())*((xcg_aft_potato -xnfree)/c_mac_w)*(1/Vtrim)
    return deriv_trim

print("------------IMPORTANT OUTPUTS FOR STABILITY-----------")
print("Stick Fixed Elevator Deflection")

print("Cmdelta_e=", Cmdelta_e(), "Control derivative is stable as Cmdelta_e < 0")
if Cmdelta_e() >0:
    print("Cmdelta_e is positive! Should be negative!!")

print((delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]))
print("Slope elevator deflection vs angle =",(delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]), "deg/deg")
if (delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]) >0:
    print("Unstability in the deflection of the elevator due to change in angle of attack :(")
else:
    print("Stability of elevator deflection due to angle of attack :)")

print("Check if close:", "Vtrim=",V_trim[100], "V_cruise=", V_cruise)
print("Slope elevator deflection vs velocity =",(delta_e(V_trim[98])-delta_e(V_trim[97]))/(V_trim[98]-V_trim[97]), "deg/m/s")
if (delta_e(V_trim[98])-delta_e(V_trim[97]))/(V_trim[98]-V_trim[97])<0:
    print("Unstability in the deflection of the elevator due to change in velocity :(")
else:
    print("Stability of elevator deflection due to velocity :)")

print("Check Control Force Stability:")
if deriv_controlForce() <0:
    print("Instability in control force :(", deriv_controlForce())
else:
    print("Stability for Control force derivative at trim speed = ", deriv_controlForce(), "N/m/s")

print("delta_te0 =", trimtab_0())
print("delta_te=", delta_te)
if delta_te > trimtab_0():
    print("Stable for delta_te0 and delta_te")
else:
    print("Unstable for the trimtab deflection")

# plt.plot(V_controlforce, delta_e(V=V_controlforce))
# plt.axvline(x=V_cruise, color="black")
# plt.grid(True)
# plt.xlabel("Velocity (m/s)")
# plt.ylabel("Elevator Deflection (degrees)")
# plt.title("Elevator Deflection Curve - Stick Fixed")
# plt.legend(["Elevator Deflection", "Vcruise"])
# plt.show()
#
# plt.plot(alpha_trim, delta_e_alpha(alpha=alpha_trim))
# plt.axvline(x=3.5, color="black")
# plt.grid(True)
# plt.xlabel("Angle of Attack (degrees)")
# plt.ylabel("Elevator Deflection (degrees)")
# plt.title("Elevator Deflection Curve - Stick Fixed")
# plt.legend(["Elevator Deflection", "alpha at cruise"])
# plt.show()
#
# plt.plot(V_controlforce, TailLoad(V = V_controlforce))
# plt.axvline(x=V_cruise, color="black")
# plt.grid(True)
# plt.xlabel("Velocity (m/s)")
# plt.ylabel("Tail Load (N)")
# plt.title("Horizontal tail load curve")
# plt.legend(["Tail loading", "Vcruise"])
# plt.show()

plt.plot(V_controlforce, ControlForce(V=V_controlforce, delta_te= delta_te))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Control Force (N)")
plt.axvline(x=V_cruise, color="black")          # V_cr
plt.axvline(x=124, color="g")                   # V_trim
plt.axvline(x=178, color="r")                   # V_max
plt.xlim(80,180)
plt.legend(["Elevator control force", "Vcruise", "Vtrim", "Vmax"])
plt.show()