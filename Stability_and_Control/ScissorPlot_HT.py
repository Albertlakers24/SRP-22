import numpy as np
import math as m
import matplotlib.pyplot as plt
from Constants.AircraftGeometry import S_w, Aw, c_mac_w,c_rw, Sweep_quarterchordw, Sweep_halfchordw, SweepLE, d_f_outer, SwfS_flap, c_accent_c_flap, b_flap, bw, taperw, l_f, cf_over_cprime_flap
from Constants.Empennage_LandingGear import lh, Vh_V, lambdahalf_h, cmac_h, Av, Sweep_halfchord_VT, Av_eff
from Constants.Aircraft_Geometry_Drawing import ln1, ln2, bn1,bn2,l_fn
from Constants.Masses_Locations import LEMAC, xcg_gear,xcg_front_potato, xcg_aft_potato
from Constants.Aerodynamics import Cm0_airfoil, CL0_Land, DeltaCLflaps, CL_DesCruise,CL_Alpha_Wing, CL_Max_Clean, CL_MaxLand, DeltaClmax
from Constants.MissionInputs import M_cruise, ISA_calculator, h_cruise, dt_cruise, V_approach, g, landing_critical, dt_land
from Constants.Masses_Locations import m_mto
from Stability_and_Control.Initial_Tail_Sizing import Geometry_HT
#TODO:  revisit xac location ? TOTAL

# Imported Variables :  Aircraft Geometry Drawings to do
v_t_w = 3.3                     # Vertical distance HT and wing     [m]
W_landing = m_mto*g             # Landing weight                    [N]

# Inputs for Geometry of HT: TO BE CHANGED IF NECESSARY
Sh = 15.0
A_h = 4.7
taperh = 0.75

# Inputs : Constants
SM = 0.05                       # Stability Margin                  [-]
eta = 0.95                      # Airfoil Efficiency Factor         [-]
kn = -0.4                       # Number                            [-]     -> nacelle mounted in front of LE
C_L_h_adj = -0.8                # CL for adjustable tail            [-]

#Supporting Equations
beta = np.sqrt(1-M_cruise**2)                           # Compressibility factor    [-]
cg = S_w/bw                                             # Mean geometric chord      [m]
Sweep_beta = (m.atan(np.tan(SweepLE)/beta)/np.pi)*180   # Sweep beta                [deg]
rho = ISA_calculator(h = h_cruise, dt=dt_cruise)[2]     # Density                   [kg/m3]

print("-----Needed for AC calculations graph (x_ac_w)-------")
print("beta*A = ", beta*Aw)
print("taper = ", taperw)
print("Sweep_beta=", Sweep_beta, "degrees")

print("-----Needed for flap contribution graphs (mu1, mu2, mu3)------")
print("Type of flap = Slotted (Fowler) flaps")
print("taper wing =", taperw)
print("b_flap/b", b_flap/(bw/2))
print("cf/c'=", cf_over_cprime_flap)

# Graph Variables       TODO: TO BE FILLED IN
x_ac_w = 0.25                   # Wing contribution to ac           [MAC]   Lecture 7 (SEAD), Slide 31 (Torenbeek)
mu1 =0.20                       # Lecture 8, Slide 20               [-]     Lecture 8 (SEAD) (Torenbeek)
mu2 =0.58                       # Lecture 8, Slide 21               [-]     Lecture 8 (SEAD) (Torenbeek)
mu3 =0.058                      # Lecture 8, Slide 21               [-]     Lecture 8 (SEAD) (Torenbeek)

def Location_in_MAC(Location):
    """
    :return: Location (MAC)
    """
    xcg_MAC = (Location-LEMAC)/c_mac_w
    return xcg_MAC

def C_L_alpha(A, lambdahalf):
    """ CHECKED
    :param A: Aspect Ratio
    :param lambdahalf: half chord sweep
    :return: C_L_alpha (per rad)
    """
    C_L_alpha = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * beta / eta)**2 * (1 + ((np.tan(lambdahalf)) ** 2 / beta ** 2)))))
    return C_L_alpha

def CL_alpha_Ah():
    """ CHECKED
    :return: CL_alpha_Ah (per rad)
    """
    Snet = S_w - c_rw*d_f_outer
    C_L_alpha_Ah = (CL_Alpha_Wing * (1 + (2.15 * (d_f_outer / bw))) * (Snet / S_w)) + ((np.pi / 2) * (d_f_outer ** 2 / S_w))
    return C_L_alpha_Ah

def Downwash():
    """
    :param mtv; phi = sin^-1(mtv/r)
    :return: Downwash gradient (-)
    """
    r = lh/(bw/2)
    mtv = v_t_w/(bw/2)
    K_EA = ((0.1124 + 0.1265 * Sweep_quarterchordw + 0.1766 * Sweep_quarterchordw ** 2) / r ** 2) + 0.1024 / r + 2
    K_EA0 = (0.1124 / r ** 2) + (0.1024 / r) + 2
    downwash = (K_EA/K_EA0)*(r/(r**2+mtv**2) * (0.4876/np.sqrt(r**2+0.6319+mtv**2)) +(1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*(CL_Alpha_Wing/(np.pi*Aw))
    # downwash = 4 / (Aw + 2)
    return downwash

print("Downwash=", Downwash())

def AC_location():
    """
    AC location (excluding the tail)
    Note: lambda_quarterchord = 0
    :param bf, hf, l_fn, c_bar (m)
    :param CL_alpha_Ah
    :param S (m^2)
    :param cg; geometric center (m)
    :param taper (-)
    :param lambda_quarterchord (radians)
    :return: x_ac (MAC)
    """
    x_ac_f1 = -((1.8 * d_f_outer * d_f_outer * l_fn) / (CL_alpha_Ah() * S_w * c_mac_w))
    x_ac_f2 = ((0.273 * d_f_outer * cg * (bw - d_f_outer)) / ((1 + taperw) * c_mac_w ** 2 * (bw + 2.15 * d_f_outer))) * np.tan(Sweep_quarterchordw)
    x_ac_n = kn*2*(bn1**2 * ln1)/(S_w * c_mac_w * CL_alpha_Ah()) + kn*2*(bn2**2 * ln2)/(S_w *c_mac_w * CL_alpha_Ah())
    x_ac = x_ac_w+x_ac_f1+x_ac_f2+x_ac_n
    #x_ac = 0.25
    return x_ac

def total_ac_calc():
    x_ac_tail = Location_in_MAC(c_mach_h*0.25 + x_h)
    eta_h = 1
    x_ac_A = AC_location()*CL_alpha_Ah() + eta_h * C_L_alpha(Ah, lambdahalf_h) * (1 - Downwash()) * (Sh / Sw) * x_ac_tail
    x_ac_tail_scaled =eta_h * C_L_alpha(Ah, lambdahalf_h) * (1 - Downwash()) * (Sh / Sw) * x_ac_tail
    x_ac_wf_scaled = AC_location()*CL_alpha_Ah()
    return x_ac_A, x_ac_tail_scaled, x_ac_wf_scaled

def flapcont(mu1, mu2, mu3):
    flap_cont_quarter = mu2 * (-mu1 * DeltaClmax * c_accent_c_flap - (CL_MaxLand + DeltaClmax * (1 - SwfS_flap)) * (1 / 8) * c_accent_c_flap*(c_accent_c_flap - 1)) + 0.7 * (Aw / (1 + 2 / Aw)) * mu3 * DeltaClmax * np.tan(Sweep_quarterchordw)
    #flap_cont_quarter = -mu1 * DeltaClmax * c_accent_c - (CL_MaxLand + DeltaClmax * (1 - SwfS)) * (1 / 8) * c_accent_c*(c_accent_c - 1)
    flap_cont = flap_cont_quarter - (CL_Max_Clean-DeltaCLflaps)*(0.25-AC_location())
    return flap_cont

def C_m_AC(mu1, mu2, mu3):
    """
    :param mu1, mu2, mu3 (-) from the graphs
    :param Cm0_airfoil (-)
    :param A; aspect ratio (-)
    :param SweepLE; (radians)
    :param bg, hf, lf, c_bar (m)
    :param CL0 (-)
    :param S (m^2)
    :param CL_alpha_Ah (rad^-1)
    :param nac_cont (cite: https://www.aerostudents.com/DSE_Reports/How_Far_Can_We_Get_Group_3_2018.pdf)
    :return:Cmac (-)
    """
    C_m_acw = Cm0_airfoil* (Aw*(np.cos(SweepLE)**2))/(Aw+2*np.cos(SweepLE))
    fus_cont = -1.8 * (1 - 2.5*d_f_outer/l_f)*(np.pi * d_f_outer * d_f_outer * l_f * CL0_Land)/(4*S_w*c_mac_w * CL_alpha_Ah())
    nac_cont = -0.05
    C_m_ac = C_m_acw + flapcont(mu1=mu1, mu2=mu2, mu3=mu3) + fus_cont + nac_cont
    #C_m_ac = -0.4
    return C_m_ac

print("Cmac=", C_m_AC(mu1=mu1, mu2=mu2, mu3=mu3))

def C_L_Ah():
    """
    CHECKED
    :param W_landing (N)
    :param rho_land: density for critical landing case (kg/m^3)
    :param S (m^2)
    :param V_landing (m/s)
    :return: CL_Ah: lift coefficien
    """
    rho_land = ISA_calculator(h=landing_critical, dt=dt_land)[2]
    C_L_Ah = 2*W_landing/(rho_land*S_w*V_approach**2)
    return C_L_Ah

# Range for the locations:
x = np.arange(-0.5,2,0.01)             #xcg/mac

#STABILITY, presented in the form y = m_s x + c_s
m_s = 1/((C_L_alpha(A = A_h, lambdahalf=lambdahalf_h) /CL_alpha_Ah())* (1-Downwash())* (lh/c_mac_w) * (Vh_V**2))
c_s = (AC_location()-0.05)*m_s
c_s_SM = AC_location()*m_s
ys = m_s * x - c_s
ys_SM = m_s*x - c_s_SM

#CONTROLLABILITY, presented in the form y = m_c x + c_c
def y_c(C_L_h):
    m_c = 1 / ((C_L_h / C_L_Ah()) * (lh / c_mac_w) * (Vh_V ** 2))
    c_c = ((C_m_AC(mu1=mu1, mu2=mu2, mu3=mu3) / C_L_Ah()) - AC_location()) * m_c
    yc = m_c * x + c_c
    return yc

# CG location range
xcg = [Location_in_MAC(xcg_front_potato), Location_in_MAC(xcg_aft_potato)]
#xcg = [Location_in_MAC(xcg_front_potato)-0.05, Location_in_MAC(xcg_aft_potato)+0.05]       # Including Margin of 0.05
ycg = [Sh/S_w, Sh/S_w]
print("------------------------------------------------------")
print("--------------- VISIBILITY CHECKS ------------------------")
print("Check Location of the landing gear is behind the cg")
print("Check the area of the horizontal tail")
print("--------------- HT Dimensions ------------------------")
print("taperh = ", taperh)
print("Sh =", Sh, "m^2")
print("Ah =", A_h)
print("bh =", Geometry_HT(Sh, A_h, taperh)[0], "m")
print("cr_h =", Geometry_HT(Sh, A_h, taperh)[1], "m")
print("ct_h =", Geometry_HT(Sh, A_h, taperh)[2], "m")

#Plotting the curves
plt.plot(x, ys_SM, 'g', linestyle = '--', label='Neutral stability Line')
plt.plot(x,ys, 'g', label="Stability Line")
plt.plot(x, y_c(C_L_h =C_L_h_adj), 'b', label = 'Controllability, for adjustable tail')
#plt.plot(x, y_c(C_L_h =C_L_h_fix), 'g', label = 'Controllability, for fixed tail')
#plt.plot(x, y_c(C_L_h =C_L_h_mov), 'r', label = 'Controllability, for full moving tail')
plt.axvline(x=Location_in_MAC(xcg_gear), color="black", label='Location Landing Gear')
plt.plot(xcg, ycg, 'r', label='CG range')
plt.xlabel('xcg [MAC]', color='#1C2833', weight="bold", fontsize=15)
plt.ylabel('Sh/S', color='#1C2833', weight="bold", fontsize=15)
plt.legend(loc='upper left', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(-0.2,1.0)
plt.ylim(0,0.6)
plt.grid(True)
plt.show()