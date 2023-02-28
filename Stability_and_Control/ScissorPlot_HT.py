import numpy as np
import math as m
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Wing_Loading_Diagram import V_approach_stall, beta_s_land_fc
from Class_I_Weight_Estimation.Class_I_fuel_cell import m_MTOW
from Initial_Aircraft_Sizing.Empennage_Design import l_h, Ah, bv, Sh, c_mach_h, x_h
from Initial_Aircraft_Sizing.Wing_planform import c_mac, M_cruise, Sw, A, b, c_r, taper
from Initial_Aircraft_Sizing.Fuselage import D_outer, l_f
from Aerodynamic_characteristics.AeroData import Cm0_AF, CL0_wing_CR, CL_Des_CR, CL_alpha
from Aerodynamic_characteristics.HLD_design_decision import DCl, c_prime_over_c_land, ld, lm, SwfS, cf_over_cprime

#TODO:  revisit xac location ? TOTAL

# Imported Variable : Aircraft Geometry
Sw = 59                         # Wing : Surface area               [m^2]
ARw = 14                        # Wing : Aspect Ratio               [-]
c_bar = 3                       # Wing : Mean aerodynamic chord     [m]
cr = 3                          # Wing : Root chord                 [m]
lambdahalf_w = -1.81*np.pi/180  # Wing : Half sweep                 [rad]
lambda_quarterchord = 0         # Wing : c/4 sweep                  [rad]
SweepLE = 1.81*np.pi/180        # Wing : Leading edge sweep         [rad]
A_h = 4                         # HT : Aspect Ratio                 [-]
lh = 12                         # HT : Tail arm                     [m]
lambdahalf_h = 2.05*np.pi/180   # HT : Half sweep                   [rad]
bf = 3                          # Fuselage : Width                  [m]
hf = bf                         # fuselage : Height                 [m]
l_fn = 10.24                    # Length nose to LE wing            [m] TODO: = LEMAC
Swf_over_S = SwfS               # Wing area by flaps                [-]
c_accent_c = c_prime_over_c_land # cprime over c                    [-]
b_flap = ld-lm                  # span flaps                        [m]

# Imported Variables :  Aircraft Geometry Drawings
v_t_w = 4.166                   # Vertical distance HT and wing     [m]
ln1 = 1.983                     # Big engine distance of front nacelle to c/4    [m]
ln2 = 1.34                      # Small engine distance of front nacelle to c/4  [m]
bn1 = 0.766                     # Big engine width of nacelle       [m]
bn2 = 0.596                     # Small engine width of nacelle     [m]

# Imported Variables : Flight Performance
M = M_cruise                    # Mach number at cruise             [-]
rho = ISA_calculator(h = h_cruise, dt=dt_cruise)[2]   # Density     [kg/m3]
V_landing = V_approach          # Velocity at landing               [m/s]    or use V_approach_stall?

# Imported Variables : Aerodynamics
Cm_0airfoil = Cm0_AF            # Cm for airfoil at alpha=0         [-]        -0.0794
CL_alpha_w = CL_alpha           # Slope lift curve wing             [rad^-1]
C_L_0 = CL0_wing_CR             # CL for alpha=0 for flapped wing   [-]        0.29971
CL = CL_Des_CR                  # CL at clean cruise condition      [-]        0.630347
C_L_h_adj = -0.8                # CL for adjustable tail            [-]
C_L_h_mov = -1                  # CL for movable tail               [-]
DeltaClmax= DCl                 # Change in airfoil Cl by flaps     [-]
CL_clean = 1.45                 # Clean configuration max CL        [-]
DeltaCLflaps = 0.94             # Change in CL due to flaps         [-]

# Imported Variables : Locations
xcg_gear = 13.4                 # Location landing gear             [m]
xcg_aft_potato = 12.6           # Aft cg location (potato plot)     [m]
xcg_front_potato = 11.7         # Front cg location (potato plot)   [m]
LEMAC =11.507                   # LE of MAC                         [m]

# Inputs : Constants
SM = 0.05                       # Stability Margin                  [-]
eta = 0.95                      # Airfoil Efficiency Factor         [-]
Vh_V = 1                        # Volume fraction for T-tail        [-]
kn = -0.4                       # Number                            [-]     -> nacelle mounted in front of LE

#Supporting Equations
beta = np.sqrt(1-M**2)                  # Compressibility factor    [-]
cg = Sw/b                               # Mean geometric chord      [m]
W_landing = m_MTOW*beta_s_land_fc*g     # Landing weight            [N]
C_L_h_fix = -0.35*A_h**(1/3)            # CL for fixed horizontal tail  [-]
Sweep_beta = (m.atan(np.tan(SweepLE)/beta)/np.pi)*180       # Sweep beta        [degrees]

print("-----Needed for AC calculations graph (x_ac_w)-------")
print("beta*A = ", beta*A)
print("taper = ", taper)
print("Sweep_beta=", Sweep_beta, "degrees")

print("-----Needed for flap contribution graphs (mu1, mu2, mu3)------")
print("Type of flap = Slotted (Fowler) flaps")
print("taper =", taper)
print("b_flap/b", b_flap/(b/2))
print("cf/c'=", cf_over_cprime)

# Graph Variables
x_ac_w = 0.25                   # Wing contribution to ac           [MAC]   Lecture 7 (SEAD), Slide 31 (Torenbeek)
mu1 =0.205                      # Lecture 8, Slide 20               [-]     Lecture 7 (SEAD) (Torenbeek)
mu2 =0.65                       # Lecture 8, Slide 21               [-]     Lecture 7 (SEAD) (Torenbeek)
mu3 =0.06                       # Lecture 8, Slide 21               [-]     Lecture 7 (SEAD) (Torenbeek)

print("-----------Check Values-----------")
def Location_in_MAC(Location):
    """
    :return: Location (MAC)
    """
    xcg_MAC = (Location-LEMAC)/c_bar
    return xcg_MAC

def C_L_alpha(A, lambdahalf):
    """ CHECKED
    :param A: Aspect Ratio
    :param lambdahalf: half chord sweep
    :return: C_L_alpha (per rad)
    """
    C_L_alpha = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * beta / eta)**2 * (1 + ((np.tan(lambdahalf)) ** 2 / beta ** 2)))))
    return C_L_alpha

print("CL_alpha_h", C_L_alpha(Ah, lambdahalf_h), "per rad")

def CL_alpha_Ah():
    """ CHECKED
    :return: CL_alpha_Ah (per rad)
    """
    Snet = Sw - cr*bf
    C_L_alpha_Ah = (CL_alpha_w * (1 + (2.15 * (bf / b))) * (Snet / S)) + ((np.pi / 2) * (bf ** 2 / S))
    return C_L_alpha_Ah

print("CL_alpha_tailless=", CL_alpha_Ah(), "per rad")

def Downwash():
    # TODO: revisit equations, in T-Tail there is probably no downwash at all.
    """
    :param mtv; phi = sin^-1(mtv/r)
    :return: Downwash gradient (-)
    """
    r = lh/(b/2)
    mtv = v_t_w
    K_EA = ((0.1124 + 0.1265 * lambda_quarterchord + 0.1766 * lambda_quarterchord ** 2) / r ** 2) + 0.1024 / r + 2
    K_EA0 = (0.1124 / r ** 2) + (0.1024 / r) + 2
    downwash = (K_EA/K_EA0)*(r/(r**2+mtv**2) * (0.4876/np.sqrt(r**2+0.6319+mtv**2)) +(1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*(CL_alpha_w/(np.pi*A_w))
    # downwash = 4 / (A + 2)
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
    x_ac_f1 = -((1.8 * bf * hf * l_fn) / (CL_alpha_Ah() * Sw * c_bar))
    x_ac_f2 = ((0.273 * bf * cg * (b - bf)) / ((1 + taper) * c_bar ** 2 * (b + 2.15 * bf))) * np.tan(lambda_quarterchord)
    x_ac_n = kn*2*(bn1**2 * ln1)/(Sw * c_bar * CL_alpha_Ah()) + kn*2*(bn2**2 * ln2)/(Sw * c_bar * CL_alpha_Ah())
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

print(C_L_alpha(Ah, lambdahalf_h), CL_alpha_Ah(), Downwash(), Sh)
print(Location_in_MAC(c_mach_h*0.25 + x_h))
print("AC location entire aircraft =", total_ac_calc(), "MAC")
print("AC location tailless aircraft= ",AC_location(), "MAC")
print("xcg_aft=", Location_in_MAC(xcg_aft_potato), "MAC")
print("xcg_front=", Location_in_MAC(xcg_front_potato), "MAC")

def flapcont(mu1, mu2, mu3):
    flap_cont_quarter = mu2 * (-mu1 * DeltaClmax * c_accent_c - (CL_max_landing + DeltaClmax * (1 - SwfS)) * (1 / 8) * c_accent_c*(c_accent_c - 1)) + 0.7 * (A / (1 + 2 / A)) * mu3 * DeltaClmax * np.tan(lambda_quarterchord)
    #flap_cont_quarter = -mu1 * DeltaClmax * c_accent_c - (CL_max_landing + DeltaClmax * (1 - SwfS)) * (1 / 8) * c_accent_c*(c_accent_c - 1)
    flap_cont = flap_cont_quarter - (CL_clean-DeltaCLflaps)*(0.25-AC_location())
    return flap_cont

print("flap cont", flapcont(mu1, mu2, mu3))
print("CL_max_land",CL_max_landing)

def C_m_AC(mu1, mu2, mu3):
    """
    :param mu1, mu2, mu3 (-) from the graphs
    :param Cm_0airfoil (-)
    :param A; aspect ratio (-)
    :param SweepLE; (radians)
    :param bg, hf, lf, c_bar (m)
    :param CL0 (-)
    :param S (m^2)
    :param CL_alpha_Ah (rad^-1)
    :param nac_cont (cite: https://www.aerostudents.com/DSE_Reports/How_Far_Can_We_Get_Group_3_2018.pdf)
    :return:Cmac (-)
    """
    C_m_acw = Cm_0airfoil* (A*(np.cos(SweepLE)**2))/(A+2*np.cos(SweepLE))
    fus_cont = -1.8 * (1 - 2.5*bf/l_f)*(np.pi * bf * hf * l_f * C_L_0)/(4*Sw*c_bar * CL_alpha_Ah())
    nac_cont = -0.05
    C_m_ac = C_m_acw + flapcont(mu1=mu1, mu2=mu2, mu3=mu3) + fus_cont + nac_cont
    #C_m_ac = -0.4
    return C_m_ac

print("Cmac=", C_m_AC(mu1=mu1, mu2=mu2, mu3=mu3))


print("cprimeoverc", c_accent_c)
print("Swf/S=", Swf_over_S)
print("Delta Clmax=", DeltaClmax)

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
    C_L_Ah = 2*W_landing/(rho_land*Sw*V_landing**2)
    return C_L_Ah
print("CL for aircraft less tail=",C_L_Ah())

x = np.arange(-0.5,2,0.01)             #xcg/mac
#Final equation for STABILITY, presented in the form y = m_s x + c_s
m_s = 1/((C_L_alpha(A = Ah, lambdahalf=lambdahalf_h) /CL_alpha_Ah())* (1-Downwash())* (lh/c_bar) * (Vh_V**2))
c_s = (AC_location()-0.05)*m_s
c_s_SM = AC_location()*m_s
print("xbarAC=",AC_location())
ys = m_s * x - c_s
ys_SM = m_s*x - c_s_SM

#CONTROLLABILITY, presented in the form y = m_c x + c_c
def y_c(C_L_h):
    m_c = 1 / ((C_L_h / C_L_Ah()) * (lh / c_bar) * (Vh_V ** 2))
    c_c = ((C_m_AC(mu1=mu1, mu2=mu2, mu3=mu3) / C_L_Ah()) - AC_location()) * m_c
    yc = m_c * x + c_c
    return yc

# CG location range
xcg = [Location_in_MAC(xcg_front_potato)-0.05, Location_in_MAC(xcg_aft_potato)+0.05]
ycg = [Sh/Sw, Sh/Sw]
print("Sh=", Sh)
#Plotting the curves
plt.plot(x, ys_SM, 'g', linestyle = '--', label='Neutral stability Line')
plt.plot(x,ys, 'g', label="Stability Line")
plt.plot(x, y_c(C_L_h =C_L_h_adj), 'b', label = 'Controllability, for adjustable tail')
#plt.plot(x, y_c(C_L_h =C_L_h_fix), 'g', label = 'Controllability, for fixed tail')
#plt.plot(x, y_c(C_L_h =C_L_h_mov), 'r', label = 'Controllability, for full moving tail')
plt.axvline(x=Location_in_MAC(xcg_gear), color="black", label='Location Landing Gear')
plt.plot(xcg, ycg, 'r', label='CG range')
# plt.title('Horizontal Tail Sizing with Longitudinal Stability and Control', fontsize=20)
plt.xlabel('xcg [MAC]', color='#1C2833', weight="bold", fontsize=15)
plt.ylabel('Sh/S', color='#1C2833', weight="bold", fontsize=15)
plt.legend(loc='upper left', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim([-0.2, 1.0])
plt.ylim([0,0.3])
plt.grid(True)
plt.show()
#x = np.linspace(-1,1,100)   #Vary this as you