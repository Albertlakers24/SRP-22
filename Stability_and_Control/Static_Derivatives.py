import numpy as np
import math as m
from Constants.AircraftGeometry import taperw, Aw, c_mac_w,delta_a_righ, delta_a_left, t_c_ratio_w, cf_over_c_aileron, y_inboard_aileron, y_outboard_aileron, bw, Sweep_quarterchordw, y_outboard_flap, y_inboard_flap, twist, dihedralw
from Constants.Aerodynamics import Cl_Alpha_WingAirfoil, CL_DesCruise, Alpha_DesCruise, CD0_CR, CL_Alpha_Wing, CL_Alpha_HT, CLH_adj, CD0_tailhCR
from Constants.MissionInputs import M_cruise
from Constants.Stability_Control import Cybeta_v, x_AC, Cnbeta, Clbeta
from Constants.Masses_Locations import xcg_front_potato
from Constants.Aircraft_Geometry_Drawing import zv
from Constants.AircraftGeometry import d_f_outer
from Constants.Empennage_LandingGear import Sweep_quarter_chord_H

# todo for YAW RATE: calculate depending on the flap we have (p.261) -> For YAW RATE calculations
delta_f = 40                                                        # [deg] flap deflection
deltacl =  0.45*0.35*np.radians(40)*cl_alpha                        # section 8.1.2.1           (p. 261)
cl_alpha = Cl_Alpha_WingAirfoil                                     # 8.1.1.2.                  (p. 247)
alpha_deltaf = deltacl / (cl_alpha * delta_f)                       # in rad/deg [?]            (p. 454)

# todo for ROLL RATE
# New names for this file:
bf = d_f_outer
CL_alpha_w_CL = CL_Alpha_Wing                                      #Wing lift curve slope at any CL    CHECKED
CL_alpha_w_CL0 =CL_Alpha_Wing                                      #wing lift curve slope at CL=0      CHECKED
CL_alpha_M_H =Cl_Alpha_WingAirfoil                                 # Airfoil lift curve at any Moment   CHECKED
CL_alpha_H_CL = CL_Alpha_HT
CL_alpha_H_CL0 = CL_Alpha_HT
CLH = CLH_adj                                                          # CLh for horizontal tail - adjustable ADSEE Lec 8 -- CHECKED
CD0H =  CD0_tailhCR                                          # CD0 for horizontal tail - from drag estimations   CHECKED

zw =


print("FILE: Static Derivatives")

"""
In this code, multiple static stability and control derivatives have been calculated
These are listed in the following order
* Aileron, delta_a                      Aileron()
* Yaw rate, r                           YawRate()
* Roll rate, p                          RollRate()
* Rudder, delta_r                       Rudder()
* Sideslip angle, beta                  Sideslip()
* Lateral Control Department Parameter  LCDP
"""

"""Calculation for the AILERON derivatives:"""

# 5 Graph Values for the Aileron Derivatives        TODO: find graph AILERON
beta_Cl_accent_delta_k = 0.6                        # Rolling moment effectiveness parameter    [rad^-1]        Graph 10.46 (p. 476)
Cl_delta_CL_deltatheory = 1                         # Correction factor for plain flap lift     [-]             Graph 8.15 (p. 262)
Cl_deltatheory = 5.2                                # Lift effectiveness of a plain flap        [rad^-1]        Graph 8.14 (p. 260)
k_a = -0.05                                         # Correlation Constant                      [-]             Graph 10.48 (p. 480)

# 6 Graph Values for the Yaw rate Derivatives       TODO: find graph YAW RATE
Delta_Clr_alpha_delta_f_outer = 0.0078              # Effect of symmetric flap deflection       [1/(rad-deg)]   Graph 10.43 (p. 463)
Delta_Clr_alpha_delta_f_inner = 0.0030              # Effect of symmetric flap deflection       [1/(rad-deg)]   Graph 10.43 (p. 463)
Cnr_CL2 = -0.02                                     # Wing Yaw Damping derivative               [rad^-1]        Graph 10.44 (p.465)
Cnr_CD0 =-0.3                                       # Drag effect                               [rad^-1]        Graph 10.45 (p.466)
Clr_CL_CL0_M0 =0.25                                 # Lifting effect                            [rad^-1]        Graph 10.41 (p. 462)
Delta_lr_epsilont = -0.017                          # Effect of wing twist                      [1/(rad-deg)]   Graph 10.42 (p. 462)

# 6 Graph Values for the Roll rate Derivatives      TODO: find graph ROLL RATE
roll_damping_parameter_w = -0.54                    # Wing roll damping parameter               [-]             Graph 10.35 (p.450)
roll_damping_parameter_H = -0.325                   # HT roll damping parameter                 [-]             Graph 10.35 (p.450)
drag_dueto_lift_roll_damping_parameter_w = -0.01    # Wing roll damping parameter D-L           [rad^-1]        Graph 10.36 (p.452)
drag_dueto_lift_roll_damping_parameter_H = -0.0125  # HT roll damping parameter D-L             [rad^-1]        Graph 10.36 (p.452)
Cnp_epsilont = 0.0075                               # Effect of wing twist                      [(rad*deg)^-1]  Graph 10.37 (p.452)
DeltaCnp_alpha_deltaf = 0.0012                      # Effect of symmetrical flap deflection     [(rad*deg)^-1]  Graph 10.38 (p.455)

# Intermediate equations AILERON
Cl_alpha_M = Cl_Alpha_WingAirfoil                   # Airfoil lift curve at any Moment          [rad^-1]
Cl_alpha_a = Cl_Alpha_WingAirfoil                   # Airfoil lift curve covered by ailerons    [rad^-1]
beta = (1-M_cruise**2)**(1/2)                       # Compressibility factor                    [-]         Eq 10.53 (p.449)
k_factor = Cl_alpha_M*(beta/(2*np.pi))              # k factor                                  [-]         Eq 10.54 (p.451)

print("------------Needed to find 5 graph coefficients for Aileron Derivatives ---------------")
print("taper =", taperw)
print("y_outboard/(b/2)= eta_o=", y_outboard_aileron/(bw/2))
print("y_inboard/(b/2)= eta_i=", y_inboard_aileron/(bw/2))
print("beta*A/k=",beta*Aw/k_factor)
print("Sweep beta=", Sweep_beta, "degrees")
print("t/c=", t_c_ratio_w)
print("cf/c for aileron=", cf_over_c_aileron)
print("Aw=", Aw)
print("Clalpha/Clalpha_theory=", Cl_Alpha_WingAirfoil/(2*np.pi))

def Aileron():      #CHECKED
    """Roskam Book:
    Cy_delta_a (p. 474)
    Cl_delta_a (p. 474)
    Cn_delta_a (p. 479)"""
    Cy_delta_a = 0

    Caccent_l_delta =(k_factor/beta)*beta_Cl_accent_delta_k                 # [rad^-1]
    Cl_delta = Cl_delta_CL_deltatheory*Cl_deltatheory                       # [rad^-1]  Airfoil
    alpha_delta_a = Cl_delta/Cl_alpha_a                                     # [rad]
    CL_delta = alpha_delta_a*Caccent_l_delta                                # [-]       Wing
    Cl_delta_a = (CL_delta/2 + CL_delta/2)*(delta_a_left - delta_a_right)   # [-]
    delta_a = 0.5*(delta_a_left+delta_a_right)                              # [rad]
    #Cl_delta_a = Cl_aileron*delta_a                                        # [old formula]

    Cn_delta_a = k_a*CL_DesCruise*Cl_delta_a
    return Cy_delta_a, Cl_delta_a, Cn_delta_a

"""Calculation for the YAW RATE derivatives:"""

# Intermediate equations YAW RATE
x_bar = AC_loc - xcg_front_potato                                   # Distance ac to cg                 [m]     todo: most critical at fron cg right?
Delta_Clr_alpha_delta_f = Delta_Clr_alpha_delta_f_outer-Delta_Clr_alpha_delta_f_inner
alpha = Alpha_DesCruise*np.pi/180                                   # Angle of attack for cruise        [rad]

print("------------Needed to find X graph coefficients for YAW RATE derivatives ---------------")
print("Aw=", Aw)
print("x_bar/MAC =", x_bar/c_mac_w)                                                     # todo: should be about 0.2
print("Sweep_quarterchord =", Sweep_quarterchordw)
print("taper=", taperw)
print("y_outboard/(b/2)=", y_outboard_flap/(bw/2))
print("y_inboard/(b/2)=", y_inboard_flap/(bw/2))
print("Cybeta_v", Cybeta_v)

def YawRate():
    """Roskam Book:
        Cy_r (p. 460)
        Cl_r (p. 460)
        Cn_r (p. 464)"""
    Cy_r = -2*Cybeta_v*(l_v*np.cos(alpha)+zv*np.sin(alpha))/b

    B = (1-(M_cruise**2)*(np.cos(Sweep_quarterchordw))**2)**(1/2)
    num = 1 + ((A*1-B**2)/(2*B*(A*B + 2*np.cos(Sweep_quarterchordw)))) + (((A*B + 2*np.cos(Sweep_quarterchordw))/(A*B + 4*np.cos(Sweep_quarterchordw)))*(np.tan(Sweep_quarterchordw))**2)
    den = 1 + (((A + 2*np.cos(Sweep_quarterchordw))/(A + 4*np.cos(Sweep_quarterchordw)))*(np.tan(Sweep_quarterchordw))**2)

    Clr_CL_CL0 = Clr_CL_CL0_M0 *num/den
    Clrw =CL_DesCruise*Clr_CL_CL0 + 0 + Delta_lr_epsilont*twist + Delta_Clr_alpha_delta_f*alpha_deltaf*delta_f
    Clrv = -(2/b**2)*(l_v*np.cos(alpha)+zv*np.sin(alpha))*(zv*np.cos(alpha)-l_v*np.sin(alpha))*Cybeta_v
    Cl_r = Clrw + Clrv

    Cnr_w = Cnr_CL2*CL_DesCruise**2 + Cnr_CD0*CD0_CR
    Cnr_v =(2/b**2)*((l_v*np.cos(alpha)+zv*np.sin(alpha))**2)*Cybeta_v
    Cn_r = Cnr_v + Cnr_w
    return Cy_r, Cl_r, Cn_r


"""Calculation for the ROLL RATE derivatives:"""

# Intermediate equations ROLL RATE
beta = (1-M_cruise**2)**(1/2)                                                               # Eq 10.53 (p.449)
k_roll_w = CL_alpha_w_CL*beta/(2*np.pi)                                                     # Eq 10.54 (p.451)
k_roll_H = CL_alpha_M_H*(beta/(2*np.pi))                                                    # Eq 10.54 (p.451)
Delta_Clp_drag_w = drag_dueto_lift_roll_damping_parameter_w*(CL_DesCruise**2)-0.125*CD0_CR  # Eq 10.56 (p.451)
Delta_Clp_drag_H = drag_dueto_lift_roll_damping_parameter_H*(CLH**2)-0.125*CD0H             # Eq 10.56 (p.451)
Dihedraleffect_w = 1-(4*zw/(bw*np.sin(dihedralw)))+12*((zw/bw)**2)*(np.sin(dihedralw))**2   # Eq 10.44 (p.451)
Dihedraleffect_HT = 1                                                                       # Eq 10.55 (p.451)

print("------------Needed to find graph coefficients for Roll Rate derivatives ---------------")
print("For wing parameters:")
print("Aw=", Aw)
print("taperw=", taperw)
print("bf/b", bf/bw)
print("for wing: beta*A/k =", (beta*Aw)/k_roll_w)
print("beta=", beta)
print("Quarter chord sweep wing =", Sweep_quarterchordw*180/np.pi, "deg")
print("for wing: Sweep_Beta =", m.atan(np.tan(Sweep_quarterchordw)/beta)*180/np.pi, "deg")

print("For HT parameters:")
print("taperH =", taperh)
print("for H-tail: beta*A/k =", (beta*Ah)/k_roll_H)
print('Quarter chord sweep Htail=', Sweep_quarter_chord_H*180/np.pi, "deg")
print("for H-tail: Sweep_Beta =", m.atan(np.tan(Sweep_quarter_chord_H)/beta)*180/np.pi, "deg")
print('Ah=', Ah)

def RollRate_Coefficient(roll_damping_parameter, k_roll, CL_alpha_i_CL, CL_alpha_i_CL0, Delta_Clp_drag_i, Dihedraleffect):
    # Eq 10.52 (p.449)
    Clpi = roll_damping_parameter * (k_roll / beta) * (CL_alpha_i_CL / CL_alpha_i_CL0)*Dihedraleffect + Delta_Clp_drag_i
    return Clpi

def RollRate():
    Cy_p = 2*Cybeta_v*(zv*np.cos(alpha)-l_v*np.sin(alpha))/bw

    Clpw = RollRate_Coefficient(roll_damping_parameter=roll_damping_parameter_w, k_roll=k_roll_w, CL_alpha_i_CL=CL_alpha_w_CL, CL_alpha_i_CL0= CL_alpha_w_CL0, Delta_Clp_drag_i=Delta_Clp_drag_w, Dihedraleffect=Dihedraleffect_w)
    Clp_h = RollRate_Coefficient(roll_damping_parameter=roll_damping_parameter_H, k_roll=k_roll_H, CL_alpha_i_CL=CL_alpha_w_CL, CL_alpha_i_CL0=CL_alpha_H_CL0, Delta_Clp_drag_i=Delta_Clp_drag_H, Dihedraleffect=Dihedraleffect_HT)
    Clph =0.5*Clp_h*(Sh/Sw)*(bh/bw)**2
    Clpv =2*Cybeta_v*(zv/bw)**2
    Cl_p = Clpw + Clph + Clpv

    Brollrate = (1-(M_cruise**2)*(np.cos(Sweep_quarterchordw))**2)**(1/2)

    part3 = A+6*(A+np.cos(Sweep_quarterchordw))*((x_bar/c_mac)*(np.tan(Sweep_quarterchordw)/A)+(np.tan(Sweep_quarterchordw)**2/12))
    part4 =A + 4*np.cos(Sweep_quarterchordw)
    Cnp_CL_CL0_M0 = -(1/6)*(part3/part4)

    part1 = (A + 4*np.cos(Sweep_quarterchordw))/(A*Brollrate + 4*np.cos(Sweep_quarterchordw))
    part2a = A*Brollrate + 0.5 * (A*Brollrate + np.cos(Sweep_quarterchordw) * np.tan(Sweep_quarterchordw) ** 2)
    part2b = A + 0.5*(A+np.cos(Sweep_quarterchordw)*np.tan(Sweep_quarterchordw)**2)
    Cnp_CL_CL0 = part1*(part2a/part2b)*Cnp_CL_CL0_M0              #Eq 10.63 (p.453) per radians
    #alpha_deltaf = deltacl/cl_alpha*delta_f

    Cnpw =Cnp_CL_CL0*CL_DesCruise + Cnp_epsilont*twist + DeltaCnp_alpha_deltaf*alpha_deltaf*delta_f
    Cnpv = -(2/b**2)*(l_v*np.cos(alpha) + zv*np.sin(alpha))*(zv*np.cos(alpha)-l_v*np.sin(alpha)-zv)*Cybeta_v
    Cn_p = Cnpw+Cnpv
    return Cy_p, Cl_p, Cn_p

def Sideslip():
    Cybeta = 1
    Clbeta = 1
    Cnbeta = 1
    return Cybeta, Clbeta, Cnbeta

def Rudder():
    Cydelta_r = Deriv_Rudder()[0]
    Cldelta_r = 1
    Cn_delta_r = Deriv_Rudder()[1]
    return Cydelta_r, Cldelta_r, Cn_delta_r

"""Calculation for the LCDP:"""
LCDP = Cnbeta-Clbeta*(Aileron()[2]/Aileron()[1])



print("-------------RESULTS STABILITY AND CONTROL DERIVATIVE--------------")
print("Aileron Design")
print("Cy_delta_a=", Aileron()[0], "Side-force-due-to-aileron derivative")
print("Cl_delta_a=", Aileron()[1], "roll control power")
print("Cn_delta_a=",Aileron()[2], "yawing moment due to aileron derivative")

print("Yaw Rate")
print("Cy_r =", YawRate()[0], "Side-force-due-to-yaw-rate derivative")
print("Cl_r =", YawRate()[1], "rolling-moment-due-to-yaw-rate derivative")
print("Cn_r = ", YawRate()[2], "yawing-moment-due-to-yaw-rate derivative")

print("Roll Rate")
print("Cy_p =", RollRate()[0], "Side-force-due-to-roll-rate derivative")
print("Cl_p =", RollRate()[1], "rolling damping derivative")
df = D_outer
if df/b < 0.3:
    print("df/b<0.3 thus Clp is normally negligible for the aircraft", df/b)
print("Cn_p = ", RollRate()[2], "yawing-moment-due-to-roll-rate derivative")

print("Sideslip Angle")
print("Cy_beta =", Sideslip()[0], "???")
print("Cl_beta =", Sideslip()[1], "?? derivative")
print("Cn_beta = ", Sideslip()[2], "?? derivative")

print("Rudder")
print("Cy_deltar =", Rudder()[0], "side-force-due-to-rudder derivative")
print("Cl_deltar =", Rudder()[1], "rolling-moment-due-to-rudder derivative")
print("Cn_deltar = ", Rudder()[2], "rudder control power")

print("Lateral Control Department Parameters")
print("LCDP=", LCDP)


# from Initial_Aircraft_Sizing.Wing_planform import M_cruise, A, taper, b, c_mac, Sw
# from Aerodynamic_characteristics.HLD_design_decision import cf_over_cprime, c_prime_over_c_to, ld, lm, DCL,Cf_C
# from Initial_Aircraft_Sizing.Empennage_Design import l_v, Sh, bh, Ah, taperh
# from Control_and_Stability.Scissorplot import AC_location, Sweep_beta
# #from Control_and_Stability.Lateral_Control_derivatives import Cybeta_v
# from Initial_Aircraft_Sizing.Fuselage import D_outer
# #from Control_and_Stability.Lateral_Control_derivatives import Cldelta_over_Cldeltatheory, Clbeta, Cn_beta
# from Aerodynamic_characteristics.AeroData import Cl_alpha_AF, CL_Des_CR, CL_alpha, Cl_alpha_AFtail, CL_alpha_htail, CL_HTail_CR
# from Aerodynamic_characteristics.Aileron_design import aileron_dA_lower, aileron_dA_upper, aileron_lm, aileron_l1, aileron_Ca_c
# from Aerodynamic_characteristics.Wing_lift_estimation import Calculate_wingsweep
# from Stability_and_Control.Vertical_Tail import Deriv_Rudder