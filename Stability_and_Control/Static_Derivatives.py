import numpy as np
import math as m
from Constants.AircraftGeometry import taperw, Aw, delta_a_righ, delta_a_left, t_c_ratio_w, cf_over_c_aileron, y_inboard_aileron, y_outboard_aileron, bw, Sweep_quarterchordw, y_outboard_flap, y_inboard_flap
from Constants.Aerodynamics import Cl_Alpha_WingAirfoil, CL_DesCruise
from Constants.MissionInputs import M_cruise
from Constants.Derivative_List import Cybeta_v

zv = 3.3026                                                         # DOROTEYA             -- CHECKED
alpha = 4 *np.pi/180                                                # [rad] cruise         -- CHECKED
delta_f = 40                                                        # [deg] flap deflection-- CHECKED
CD0w =0.008216348148174173                                          # Imported from MEGHA  -- CHECKED
xcg = 12.31402435                                                   # Most aft CG - from nose of aircraft  -- CHECKED
AC_loc = (0.35*c_mac)  + 11.507   #AC_location()                    # TODO: Get correct AC location - front of A/C
x_bar = xcg- AC_loc
Delta_Clr_alpha_delta_f = Delta_Clr_alpha_delta_f_outer-Delta_Clr_alpha_delta_f_inner   # -- CHECKED
epsilon_t =0                                                        #No twist --  CHECKED
dihedral = 0                                                        # -- CHECKED
cl_alpha =  Cl_alpha_M * np.sqrt(1-M_cruise**2)                     #lift-curve slope airfoil -- CHECKED
deltacl =  0.45*0.35*np.radians(40)*cl_alpha                        # section 8.1.2.1 -- CHECKED
wing_sweep_quarter = 0                                              # CHECKED
#Cl_alpha_a = Cl_alpha_M
alpha_deltaf = deltacl / (cl_alpha * delta_f)                       # in rad/deg -- CHECKED




print("FILE: Static Derivatives")

"""
In this code, multiple static stability and control derivatives have been calculated
These are listed in the following order
* Aileron, delta_a                      Aileron()
* Yaw rate, r                           YawRate()
* Rudder, delta_r                       Rudder()
* Roll rate, p                          RollRate()
* Lateral Control Department Parameter  LCDP
* Sideslip angle, beta                  Sideslip()
"""

"""Calculation for the AILERON derivatives:"""

# 5 Graph Values for the Aileron Derivatives        TODO: find graph AILERON
beta_Cl_accent_delta_k = 0.6                        # Rolling moment effectiveness parameter    [rad^-1]    Graph 10.46 (p. 476)
Cl_delta_CL_deltatheory = 1                         # Correction factor for plain flap lift     [-]         Graph 8.15 (p. 262)
Cl_deltatheory = 5.2                                # Lift effectiveness of a plain flap        [rad^-1]    Graph 8.14 (p. 260)
k_a = -0.05                                         # Correlation Constant                      [-]         Graph 10.48 (p. 480)

# 6 Graph Values for the Yaw rate Derivatives       TODO: find graph YAW RATE
Delta_Clr_alpha_delta_f_outer = 0.0078              # Effect of symmetric flap deflection       [1/(rad-deg)]   Graph 10.43 (p. 463)
Delta_Clr_alpha_delta_f_inner = 0.0030              # Effect of symmetric flap deflection       [1/(rad-deg)]   Graph 10.43 (p. 463)
Cnr_CL2 = -0.02                                     # Wing Yaw Damping derivative               [rad^-1]        Graph 10.44 (p.465)
Cnr_CD0 =-0.3                                       # Drag effect                               [rad^-1]        Graph 10.45 (p.466)
Clr_CL_CL0_M0 =0.25                                 # Lifting effect                            [rad^-1]        Graph 10.41 (p. 462)
Delta_lr_epsilont = -0.017                          # Effect of wing twist                      [1/(rad-deg)]   Graph 10.42 (p. 462)

# Intermediate equations AILERON
Cl_alpha_M = Cl_Alpha_WingAirfoil                   # Airfoil lift curve at any Moment          [rad^-1]
Cl_alpha_a = Cl_Alpha_WingAirfoil                   # Airfoil lift curve covered by ailerons    [rad^-1]
beta = (1-M_cruise**2)**(1/2)                       # Compressibility factor                    [-]
k_factor = Cl_alpha_M*(beta/(2*np.pi))              # k factor                                  [-]         eq. 10.54 (p. 451)

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

print("------------Needed to find X graph coefficients for YAW RATE derivatives ---------------")
print("Aw=", Aw)
print("x_bar/MAC =", x_bar/c_mac)       # todo: should be about 0.2
print("Sweep_quarterchord =", Sweep_quarterchordw)
print("taper=", taperw)
print("y_outboard/(b/2)=", y_outboard_flap/(bw/2))
print("y_inboard/(b/2)=", y_inboard_flap/(bw/2))
print("Cybeta_v", Cybeta_v)
print('alpha delta f', alpha_deltaf)

def YawRate():
    """Roskam Book:
        Cy_r (p. )
        Cl_r (p. )
        Cn_r (p. )"""
    Cy_r = -2*Cybeta_v*(l_v*np.cos(alpha)+zv*np.sin(alpha))/b

    B = (1-(M_cruise**2)*(np.cos(Sweep_quarterchordw))**2)**(1/2)
    num = 1 + ((A*1-B**2)/(2*B*(A*B + 2*np.cos(wing_sweep_quarter)))) + (((A*B + 2*np.cos(wing_sweep_quarter))/(A*B + 4*np.cos(wing_sweep_quarter)))*(np.tan(wing_sweep_quarter))**2)
    den = 1 + (((A + 2*np.cos(wing_sweep_quarter))/(A + 4*np.cos(wing_sweep_quarter)))*(np.tan(wing_sweep_quarter))**2)

    Clr_CL_CL0 = Clr_CL_CL0_M0 *num/den
    Clrw =CL_DesCruise*Clr_CL_CL0 + 0 + Delta_lr_epsilont*epsilon_t + Delta_Clr_alpha_delta_f*alpha_deltaf*delta_f
    Clrv = -(2/b**2)*(l_v*np.cos(alpha)+zv*np.sin(alpha))*(zv*np.cos(alpha)-l_v*np.sin(alpha))*Cybeta_v
    Cl_r = Clrw + Clrv

    Cnr_w = Cnr_CL2*CL_DesCruise**2 + Cnr_CD0*CD0w
    Cnr_v =(2/b**2)*((l_v*np.cos(alpha)+zv*np.sin(alpha))**2)*Cybeta_v
    Cn_r = Cnr_v + Cnr_w
    return Cy_r, Cl_r, Cn_r


#Imported Variables for Roll rate derivatives
bf = D_outer

print("------------Needed to find graph coefficients for Roll Rate derivatives ---------------")
print("Aw=", Aw)
print("taper=", taperw)
Cnp_epsilont =  0.0075                                              #from graph 10.37 (p. 452) -> not important if twist=0
print("bf/b", bf/b, D_outer)                                        #CHECKED
DeltaCnp_alpha_deltaf =  0.0012                                     #(rad*deg)^-1 from graph 10.38 (p. 423)  CHECKED

##Wing Parameters
print("for wing: beta*A/k =", (beta*A)/k_factor)
print("beta=", beta)
print("A=", A)
print("k =", k_factor)
print("for wing: Sweep_Beta =", m.atan(np.tan(Sweep_quarterchordw)/beta))
CL_alpha_w_CL = CL_alpha                                            #Wing lift curve slope at any CL    CHECKED
CL_alpha_w_CL0 =CL_alpha                                            #wing lift curve slope at CL=0      CHECKED
roll_damping_parameter_w = -0.54                                    #graph 10.35 (p.418)                CHECKED
k_roll_w = CL_alpha_w_CL*beta/(2*np.pi)                             # CHECKED
print('K roll factor for wing', k_roll_w)
drag_dueto_lift_roll_damping_parameter_w = -0.01                    #graph 10.36 (p. 420)   CHECKED
Delta_Clp_drag_w = drag_dueto_lift_roll_damping_parameter_w*(CL_DesCruise**2)-0.125*CD0w     #CHECKED

##Horizontal Tail Parameters
CL_alpha_M_H =Cl_alpha_AFtail                                       #Airfoil lift curve at any Moment   CHECKED
k_roll_H = CL_alpha_M_H*(beta/(2*np.pi))                            #eq. 10.54 (p. 451) CHECKED
print("taperH =", taperh)
print("for H-tail: beta*A/k =", (beta*Ah)/k_roll_H)
Sweep_quarter_chord_H = Calculate_wingsweep(6.14, 0,  0.25, taperh) # LE Sweep from Gabriel CHECKED
print('Quarter chord sweep Htail', Sweep_quarter_chord_H)
print("for H-tail: Sweep_Beta =", m.atan(np.tan(Sweep_quarter_chord_H)/beta))
print('K roll factor for H Tail', k_roll_H)
print('Ah', Ah)
roll_damping_parameter_H = -0.325                                   # graph 10.35 (p.450)   CHECKED
CL_alpha_H_CL = CL_alpha_htail                                      # Horizontal Tail lift curve slope at any CL    CHECKED
CL_alpha_H_CL0 = CL_alpha_htail                                     # Horizontal Tail lift curve slope at CL=0       CHECKED
CLH = -0.8                                                          # CLh for horizontal tail - adjustable ADSEE Lec 8 -- CHECKED
CD0H =  0.0741993961878728                                          # CD0 for horizontal tail - from drag estimations   CHECKED
drag_dueto_lift_roll_damping_parameter_H = -0.0125                  # graph 10.36 (p. 452) CHECKED
Delta_Clp_drag_H = drag_dueto_lift_roll_damping_parameter_H*(CLH**2)-0.125*CD0H
Clp_CL0_dihedral0 = roll_damping_parameter_w*k_roll_w/beta

##Imported Values for LCDP
Cnbeta = 0.09237746426846452                                        #DOROTEYA
Clbeta = -0.031725211066274815                                                          #DOROTEYA

LCDP = Cnbeta-Clbeta*(Aileron()[2]/Aileron()[1])

def RollRate_Coefficient(roll_damping_parameter, k_roll, CL_alpha_i_CL, CL_alpha_i_CL0, Delta_Clp_drag_i):
    #Assumption: no dihedral
    Clpi = roll_damping_parameter * (k_roll / beta) * (CL_alpha_i_CL / CL_alpha_i_CL0) + Delta_Clp_drag_i
    return Clpi

def RollRate():             # CHECKED
    Cy_p = 2*Cybeta_v*(zv*np.cos(alpha)-l_v*np.sin(alpha))/b #+ 3*np.sin(dihedral)*(1 - (4*zv/b)*np.sin(dihedral))*Clp_CL0_dihedral0

    Clpw = RollRate_Coefficient(roll_damping_parameter=roll_damping_parameter_w, k_roll=k_roll_w, CL_alpha_i_CL=CL_alpha_w_CL, CL_alpha_i_CL0= CL_alpha_w_CL0, Delta_Clp_drag_i=Delta_Clp_drag_w)
    Clp_h = RollRate_Coefficient(roll_damping_parameter=roll_damping_parameter_H, k_roll=k_roll_H, CL_alpha_i_CL=CL_alpha_w_CL, CL_alpha_i_CL0=CL_alpha_H_CL0, Delta_Clp_drag_i=Delta_Clp_drag_H)
    Clph =0.5*Clp_h*(Sh/Sw)*(bh/b)**2
    Clpv =2*Cybeta_v*(zv/b)**2
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

    Cnpw =Cnp_CL_CL0*CL_DesCruise + Cnp_epsilont*epsilon_t + DeltaCnp_alpha_deltaf*alpha_deltaf*delta_f
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



print("-------------RESULTS STABILITY AND CONTROL DERIVATIVE--------------")
print("Aileron Design")
print("Cy_delta_a=", Aileron()[0], "Side-force-due-to-aileron derivative")
print("Cl_delta_a=", Aileron()[1], "roll control power")
print("Cn_delta_a=",Aileron()[2], "yawing moment due to aileron derivative")

print("Yaw Rate")
print("Cy_r =", YawRate()[0], "Side-force-due-to-yaw-rate derivative")
print("Cl_r =", YawRate()[1], "rolling-moment-due-to-yaw-rate derivative")
print("Cn_r = ", YawRate()[2], "yawing-moment-due-to-yaw-rate derivative")

print("Lateral Control Department Parameters")
print("LCDP=", LCDP)

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