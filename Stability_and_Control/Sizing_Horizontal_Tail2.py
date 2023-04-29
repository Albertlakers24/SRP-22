import numpy as np
import matplotlib.pyplot as plt
import math as m
from Constants.MissionInputs import M_cruise, V_approach, V_cruise, g, ISA_calculator, h_cruise,dt_cruise
from Constants.Masses_Locations import m_mto, xcg_aft_potato, LEMAC
from Constants.AircraftGeometry import S_w, c_mac_w
from Constants.Empennage_LandingGear import lh,Sh,Vh_V, A_h, Sweep_quarter_chord_HT, bh, t_over_c_HT, c_th,c_rh
from Constants.Aerodynamics import Cl_Alpha_HT_Airfoil, R_HT, CL_Alpha_Wing, downwash, CL_Alpha_HT, Cmac

# To be imported from somewhere todo to be determined
V_max = 178
AlphaCL0_CR = 2*np.pi/180                   # Angle of attack at CL0        [rad]
i_h = -2*np.pi/180                          # Indince angle                 [rad]

# Values for this file
CNh_delta = 0.024*180/np.pi                 # Normal force gradient         [rad^-1]    Graph (pg. 317 FD reader)
x_ac_w_MAC = 0.25                           # Location ac of the main wing  [MAC]

# DESIGN CHOICES ELEVATOR
delta_te = -11*np.pi/180                    # For stable stick free aircraft: delta_te > delta_te0          [rad]            Initial: -11 deg
deriv_deltae_se = 2.2                       #                               [rad/m]         normal values: 2.2-26 rad/m (1.25-1.5deg/cm) --- pg 400 FD reader
inboard_elevator = 0.3                      # Inboard location elevator     [m]
outboard_elevator = bh/2                    # Outboard location elevator    [m]
cb_over_cf = 0.3                            # Hinge line position elevator  [-]
cf_over_c = 0.3                             # Relative elevator chord       [-]
delta_f_elevator = 30                       # Elevator flap deflection      [deg]
Sweep_hl = 11*np.pi/180                     # Sweep at the hingline         [rad] todo determine with drawing
tc_over_two_cf = 0.2485                     # tc = thickness chord at hingeline and cf = chord of flap (p. 509)

# SIZING ELEVATOR
angle = m.atan((c_rh-c_th)/(bh/2))                      # Angle                                 [rad]
root_H = c_rh-inboard_elevator*np.tan(angle)            # Chord of H at start elevator          [m]
y = np.tan(angle)*(bh/2 - outboard_elevator)            # Difference change                     [m]
Ct_e = (c_th+y)*cf_over_c                                                                               # Tip chord elevator        [m]
Cr_e = root_H*cf_over_c                                                                                 # Root chord elevator       [m]
MACe = Cr_e*0.75                                                                                        # MAC elevator              [m] todo estimated
Se = Ct_e*(outboard_elevator-inboard_elevator)+0.5*(Cr_e-Ct_e)*(outboard_elevator-inboard_elevator)     # Elevator surface area     [m^2]
Eta_i= inboard_elevator/(bh/2)
Eta_o= outboard_elevator/(bh/2)

""""--------------------------DETERMINE GRAPH VALUES------------------------------------"""
print("Needed for Hingemoment_Coefficients 3D calculations for the ELEVATOR")
print("A_h =", A_h)
print("Eta_i =", inboard_elevator/(bh/2))
print("Eta_o =", outboard_elevator/(bh/2))
print("cb'/cf'=", cb_over_cf*cf_over_c)
print("cf'/c'=", cf_over_c)
print("delta_f =", delta_f_elevator, "deg")
print("cf/c =", cf_over_c)
print("t/c_HT =", t_over_c_HT)

# Determined values from Graphs for Hingemoment_Coefficients 3D:
DeltaCha_over_clalphaBK = 0.01              # Graph 10.77a (p. 515)
Kalpha_i = 1.1                              # Graph 10.77b (p. 515)
Kalpha_o = 4.2                              # Graph 10.77b (p. 515)
B2 = 1.25                                   # Graph 10.77c (p. 515)
DeltaChd_cldBKd = 0.016                     # Graph 10.78a (p. 517)
Kdelta_i = 1.05                             # Graph 10.78b (p. 517)
Kdelta_o = 4.2                              # Graph 10.78b (p. 517)
alpha_d = 0.51                              # Graph 8.17   (p. 262)
cl_delta = 5.4                              # Graph 8.14   (p. 260)

# Determined values from Graphs for Hingemoment_Coefficients Airfoil:
cla_over_cla_theory_HT = 0.76               # Graph 10.64a (p. 501)
Ch_alpha_bal_over_Chalpha = 0.85            # Graph 10.65a (p. 503)
c_accent_alpha_over_chalpha_theory_HT = 0.3 # Graph 10.63a (p. 500)
chalpha_theory = -0.48                      # Graph 10.63b (p. 500)
c_accent_h_over_ch_delta_theory = 0.75      # Graph 10.69a (p. 507)
ch_delta_theory = -0.88                     # Graph 10.69b (p. 507)
Ch_delta_bal_over_Chdelta = 0.75            # Graph 10.71  (p. 509)

print("Needed for Airfoil hingemoment calculations for the ELEVATOR")
print("tan(0.5*phi_te)=", t_over_c_HT, "=t/c")
print("Reynolds Number for HT=", R_HT)
print("Nose Type = Round Nose")
print("Balance Ratio =", np.sqrt((cb_over_cf)**2-(tc_over_two_cf)**2))
print("cla_over_cla_theory_HT=",cla_over_cla_theory_HT)
print("cf/c =", cf_over_c)
print("t/c =", t_over_c_HT)

cf_c_tab = 0.1                              # ct/cf, compared to flap

print("Needed for Airfoil hingemoment claculation of the TRIMTAB")
print("cf/c_tab = ", cf_c_tab)
print("cf/c =", cf_over_c)

chdeltat = -0.009                           # Graph 10.72   (p.511)
ch_cl_t = -0.058                            # Graph 10.73   (p.511)
cla_dt = Cl_Alpha_HT_Airfoil                # 8.1.1.2
alpha_dt = -0.6                             # Graph 10.74   (p.512)

def HingemomentAF_Coefficients():
    """ TRIMTAB
    :return:
    Ch_alpha_AF (p.499)
    Ch_delta_AF (p.506)
    Ch_delta_t_AF (p.510)
    """
    c_accent_h_alpha = c_accent_alpha_over_chalpha_theory_HT * chalpha_theory
    Ch_alpha_bal = c_accent_h_alpha * Ch_alpha_bal_over_Chalpha
    Ch_alpha_AF = Ch_alpha_bal/((1-M_cruise**2)**0.5)

    c_accent_h_delta = c_accent_h_over_ch_delta_theory*ch_delta_theory
    Ch_delta_bal = c_accent_h_delta*Ch_delta_bal_over_Chdelta
    Ch_delta_AF =Ch_delta_bal/((1-M_cruise**2)**0.5)

    Ch_delta_t_AF= chdeltat-ch_cl_t*cla_dt*alpha_dt

    return Ch_alpha_AF, Ch_delta_AF, Ch_delta_t_AF

def Hingemoment_Coefficients_elevator():
    """
    Coefficients due to the hingemoment for the elevator 3D:
    Ch_alpha (p.513)
    Ch_delta (p.516)
    Ch_delta can be reduced to the correct horn design. The number 0.26 is based on literature research from NASA
    Ch_delta_t (p.517) -> previously calculated
    """

    ch_alpha_M = HingemomentAF_Coefficients()[0]
    ch_delta_M = HingemomentAF_Coefficients()[1]
    ch_delta_t_M = HingemomentAF_Coefficients()[2]

    Kalpha = Kalpha_i*(1-Eta_i)-Kalpha_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaC_h_alpha = DeltaCha_over_clalphaBK*(Cl_Alpha_HT_Airfoil*B2*Kalpha*np.cos(Sweep_quarter_chord_HT))
    Ch_alpha = ((A_h*np.cos(Sweep_quarter_chord_HT))/(A_h+2*np.cos(Sweep_quarter_chord_HT)))*ch_alpha_M + DeltaC_h_alpha

    Kdelta = Kdelta_i*(1-Eta_i)-Kdelta_o*((1-Eta_o)/(Eta_o-Eta_i))
    DeltaCh_delta =DeltaChd_cldBKd*(cl_delta*B2*Kdelta*np.cos(Sweep_quarter_chord_HT)*np.cos(Sweep_hl))
    Ch_delta = np.cos(Sweep_quarter_chord_HT)*np.cos(Sweep_hl)*(ch_delta_M +alpha_d*ch_alpha_M)*((2*np.cos(Sweep_quarter_chord_HT))/(A_h+2*np.cos(Sweep_quarter_chord_HT)))+DeltaCh_delta
    Ch_delta_horn = Ch_delta*0.2

    Ch_delta_t = np.cos(Sweep_quarter_chord_HT)*np.cos(Sweep_hl)*(ch_delta_t_M+alpha_dt*ch_alpha_M*(2*np.cos(Sweep_quarter_chord_HT)/(A_h+2*np.cos(Sweep_quarter_chord_HT))))+DeltaCh_delta
    # Ch_delta_t = -0.0906472545474966

    return Ch_delta_horn, Ch_alpha, Ch_delta, Ch_delta_t

def Cmdelta_e():
    """
    Control derivative due to elevator deflection equations from FD reader p. 196
    :param CNh_delta: (rad^-1)
    :return: Cmdelta_e: (rad^-1)

    Normal value around -1.0 to -1.5
    """
    Cmdelta_e = -CNh_delta*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))
    return Cmdelta_e

def trimtab_0():
    """ CHECKED
    :return: delta_te0  (rad)
    """

    Ch_delta = Hingemoment_Coefficients_elevator()[2]
    Ch_delta_t = Hingemoment_Coefficients_elevator()[3]
    delta_te0 = ((Ch_delta /Cmdelta)*Cmac) + ((Ch_delta/CNh_delta)*CNhalpha_free*(AlphaCL0_CR+i_h)) / Ch_delta_t
    return delta_te0

# SUPPORTING EQUATIONS
MTOW = m_mto * g
rho = ISA_calculator(h=h_cruise, dt=dt_cruise)[2]

# Location neutral point stick free eq. 7.7 Sam1    [m]         normal value 0.483*MAC
x_ac_w_meters = c_mac_w*x_ac_w_MAC+LEMAC                                                                        # Location aerodynamic center wing                  [m]
Chalpha_Chdelta = Hingemoment_Coefficients_elevator()[1]/Hingemoment_Coefficients_elevator()[0]
CNhalpha_free= CL_Alpha_HT - CNh_delta*Chalpha_Chdelta                                                              # Normal force gradient eq. 7.5 Sam1                [rad^-1]
xnfree = ((CNhalpha_free/CL_Alpha_Wing)*(1-downwash)*(Vh_V**2)*((Sh*lh)/(S_w*c_mac_w))*c_mac_w) + x_ac_w_meters
Cmdelta = -CNh_delta*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))                                                            # eq. 5.21 Sam1                                     [rad^-1]

def ControlForce_HT(V,delta_te):
    """ For Cruise Conditions
    Control Force should be between the following values: XX < Fe < XX
    Control Curve, used for elevator control force stability (dFe/dV)Fe=0 >0
    :param: Se; Use surface area of 1 elevator (m^2) - same when you have St
    :return: Tail Load (N)
    """
    Ch_delta = Hingemoment_Coefficients_elevator()[0]
    Chdelta_t = Hingemoment_Coefficients_elevator()[3]          # To be written still!

    F_velocity_independent = (MTOW/S_w)*(Ch_delta/Cmdelta_e())*((xcg_aft_potato -xnfree)/c_mac_w)
    F_velocity_dependent= 0.5*rho*V**2*Chdelta_t*(delta_te-trimtab_0())
    a = deriv_deltae_se*Se*MACe*(Vh_V)**2
    Fe = a*(F_velocity_independent - F_velocity_dependent)
    Fe_ex_trimtab = a*F_velocity_independent
    return Fe, Fe_ex_trimtab

def ControlForceGraph(delta_te):
    V_controlforce = np.arange(0,250,0.1)
    plt.plot(V_controlforce, ControlForce_HT(V = V_controlforce,delta_te=delta_te)[0])
    plt.grid(True)
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Control Force (N)")
    plt.axvline(x=V_cruise, color="black")
    plt.axvline(x=V_approach, color="g")
    plt.axvline(x=V_max, color="r")
    plt.axhline(y=35*g, color="r")
    plt.axhline(y=-35*g,color="r")
    plt.axvspan(V_max,V_max+30,alpha=0.2, color="r")
    plt.axhspan(-35*g, -400, alpha=0.2, color="r")
    plt.axhspan(35*g,400, alpha=0.2,color="r")
    plt.xlim(0, V_max+30)
    plt.ylim(-400,100)
    plt.legend(["Elevator control force", "Vcruise", "V_approach", "Limits"])
    plt.show()
    return

ControlForceGraph(delta_te)


print("------------------------ Printed Derivatives ------------------------")
print("Chdelta =", Hingemoment_Coefficients_elevator()[0])
print("Chdelta (before horn) = ", Hingemoment_Coefficients_elevator()[2])
print("Chalpha =", Hingemoment_Coefficients_elevator()[1])
print("Chdelta_t =", Hingemoment_Coefficients_elevator()[3])
print("Cnhalpha_free =", CNhalpha_free)
print("Cmdelta =", Cmdelta)

print("--------- Elevator Requirements ----------")
if Cmdelta < 0:
    print("Elevator efficiency requirement met:")
    print("Cmdelta =", Cmdelta, "<0")
else:
    print("Elevator efficiency NOT met: Cmdelta = ", Cmdelta, "<0")

if Hingemoment_Coefficients_elevator()[0] < 0:
    print("Hinge moment derivative met:")
    print("Chdelta =", Hingemoment_Coefficients_elevator()[0],"<0")
else:
    print("Hinge moment derivative NOT met:")
    print("Chdelta =", Hingemoment_Coefficients_elevator()[0], ">0")