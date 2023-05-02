import numpy as np
import math as m
from Constants.MissionInputs import rho_0, g, V_cruise, phi, Psi
from Constants.Masses_Locations import m_mto, xcg_aft_potato, LEMAC
from Constants.AircraftGeometry import S_w, Aw, c_mac_w, bw
from Constants.Aerodynamics import CL_DesCruise, CL_Alpha_Wing, downwash, CL_Alpha_HT
from Constants.Empennage_LandingGear import Vh_V, Sh, lh
from Constants.Stability_Control import CNh_delta, Chdelta, Chalpha, CNhalpha_free

print("FILE: Derivative_Calculations_Dynamic")
gamma0 =0                               # Steady horizontal flight FD p163

Oswald = 1/((np.pi)*Aw*Psi+(1/phi))      #-      Oswald Efficiency Factor

def Zero_Derivatives():
    """Derivatives in Steady Flight
    :param W; weight [N]
    :param theta0; [rad]"""
    CX0 = m_mto*g * np.sin(gamma0) / (0.5 * rho_0 * V_cruise ** 2 * S_w)
    CZ0 =-m_mto*g * np.cos(gamma0) / (0.5 * rho_0 * V_cruise ** 2 * S_w)
    Cm0 = 0
    return CX0, CZ0, Cm0

def Velocity_Derivatives(): #CHECKED
    CXu = 2*CL_DesCruise* np.tan(gamma0)
    CZu = -2* CL_DesCruise
    Cmu = 0
    return CXu, CZu, Cmu

def Attack_Derivative(): #CHECKED
    CXalpha = CL_DesCruise*(1-((2*CL_Alpha_Wing)/(np.pi*Aw*Oswald)))
    CZalpha = -CL_Alpha_Wing
    # Cmalpha = (CL_Alpha_Wing *(x_cg - x_w)/c_mac) - CL_Alpha_HT*(1-downwash)*Vh_V**2*Sh*l_h/Sw/c_mac
    # Cmalpha1 = CL_Alpha_Wing * (x_cg - x_w)/c_mac + -((CL_Alpha_HT *(1-downwash))+(0.024*180/np.pi)*-0.8741417)*((Sh*l_h)/(Sw*c_mac))
    x_w = LEMAC + 0.25 * c_mac_w
    Cmalpha = CL_Alpha_Wing * (xcg_aft_potato - x_w)/c_mac_w - (CNhalpha_free*(1-downwash)*Vh_V**2*(Sh*lh)/(S_w*c_mac_w))

    CZalpha_dot = -CL_Alpha_HT*(Vh_V**2)*downwash*(Sh*lh/(S_w*c_mac_w))
    Cmalpha_dot = -CL_Alpha_HT * Vh_V**2 * downwash * Sh * (lh ** 2 / (S_w * c_mac_w**2))
    CXalpha_dot = 0                                                                         #unknown
    return CXalpha, CZalpha, Cmalpha, CZalpha_dot, Cmalpha_dot, CXalpha_dot #, Cmalpha1, Cmalpha2


def Pitch_Derivative(): #CHECKED
    CXq =0
    CZq = -2*CL_Alpha_HT*(Vh_V**2)*(Sh*lh/(S_w*c_mac_w))
    Cmq = -1.1 * CL_Alpha_HT * Vh_V ** 2 * Sh * lh**2 / (S_w * c_mac_w**2)
    return CXq, CZq,Cmq

def Delta_e_Derivatives(): #CHECKED
    CXde = 0
    CZde = - CNh_delta *Vh_V**2*Sh/S_w
    Cmde = -1.11
    return CXde, CZde, Cmde

xcg_aircraft = xcg_aft_potato       # m
zcg_aircraft = 3.302                # m

def Iyy(mass,xcg,zcg):
    Iyy = mass*((xcg-xcg_aircraft)**2+(zcg-zcg_aircraft))
    return Iyy

total_Iyy = 1

Ky_total = np.sqrt(total_Iyy/m_mto)
KY = Ky_total/2.34


print("--------Important Print Statements---------")

print("Oswald factor=",Oswald)
print("CX0    =", Zero_Derivatives()[0])
print("CXu    =", Velocity_Derivatives()[0])
print("CXa    =", Attack_Derivative()[0] )                                            # Positive   [FD lecture notes]
print("CXadot =", Attack_Derivative()[5])
print("CXq    =",  Pitch_Derivative()[0])
#print("CXde   =", -0.03728, "TBD")

print("CZ0    =", Zero_Derivatives()[1])
print("CZu    =", Velocity_Derivatives()[1])
print("CZa    =", Attack_Derivative()[1])
print("CZadot =", Attack_Derivative()[3])
print("CZq    =", Pitch_Derivative()[1])
#print("CZde   =", -0.69612, "TBD")

print("Cmu    =", Velocity_Derivatives()[2])
print("Cmadot =", Attack_Derivative()[4])
print("Cmq    =", Pitch_Derivative()[2])
print("Cm alpha =", Attack_Derivative()[2])

print("muc =", m_mto/(rho_0 * S_w * c_mac_w))
print("mub = ",m_mto/(rho_0 * S_w * bw))

KY2 = KY**2
KX2 = 0.002
KZ2 = 0.042
# print("KX2 = ",0.019 , "TBD")                          # squared Non-dimensional radius of gyration about the X-axis   [-]
# print("KY2 = ",1.25 * 1.114, "TBD")                    # squared Non-dimensional radius of gyration about the Y-axis   [-]
# print("KZ2 = ",0.042 , "TBD")                          # squared Non-dimensional radius of gyration about the Z-axis   [-]
# print("KXZ = ",0.002  , "TBD")                         # Non-dimensional product of inertia                            [-]

print("CXde =", Delta_e_Derivatives()[0])
print("CZde =", Delta_e_Derivatives()[1])
print("Cmde =", Delta_e_Derivatives()[2])