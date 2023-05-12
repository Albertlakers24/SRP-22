import numpy as np
import math as m
from Constants.MissionInputs import rho_0, g, V_cruise, phi, Psi,rho_cruise
from Constants.Masses_Locations import m_mto, xcg_aft_potato, LEMAC, m_oem
from Constants.AircraftGeometry import S_w, Aw, c_mac_w, bw
from Constants.Aerodynamics import CL_DesCruise, CL_Alpha_Wing, downwash, CL_Alpha_HT
from Constants.Empennage_LandingGear import Vh_V, Sh, lh
from Constants.Stability_Control import CNh_delta, Chdelta, Chalpha, CNhalpha_free

print("FILE: Derivative_Calculations_Dynamic")
gamma0 =0                               # Steady horizontal flight FD p163

Oswald = 1/((np.pi)*Aw*Psi+(1/phi))     #-      Oswald Efficiency Factor
DeltaOswald_TO = 0.0026*15              #-      Change in Oswald Efficiency Factor at TO, with 15degrees flap deflection
print("Oswald =", Oswald)
print("Oswald_TO=", DeltaOswald_TO+Oswald)

def Zero_Derivatives():
    """Derivatives in Steady Flight
    :param W; weight [N]
    :param theta0; [rad]"""
    CX0 = m_mto*g * np.sin(gamma0) / (0.5 * rho_cruise * V_cruise ** 2 * S_w)
    CZ0 =-m_mto*g * np.cos(gamma0) / (0.5 * rho_cruise * V_cruise ** 2 * S_w)
    Cm0 = 0
    return CX0, CZ0, Cm0

def Velocity_Derivatives():
    CXu = 2*CL_DesCruise* np.tan(gamma0)
    CZu = -2* CL_DesCruise
    Cmu = 0
    return CXu, CZu, Cmu

def Attack_Derivative(): #CHECKED
    CXalpha = CL_DesCruise*(1-((2*CL_Alpha_Wing)/(np.pi*Aw*Oswald)))
    CZalpha = -CL_Alpha_Wing-CNhalpha_free*(1-downwash)*Vh_V**2*(Sh/S_w)
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
    Cmq = -1.4* CL_Alpha_HT * Vh_V ** 2 * Sh * lh**2 / (S_w * c_mac_w**2)
    return CXq, CZq,Cmq
print("Cmq=", Pitch_Derivative()[2])
def Delta_e_Derivatives(): #CHECKED
    CXde = 0
    CZde = - CNh_delta *Vh_V**2*Sh/S_w
    Cmde = -1.11
    return CXde, CZde, Cmde

Xcg = xcg_aft_potato       # [m]
Ycg = 0                    # [m]

"Fuselage cg locations"
ycgfus = 0 * 0.001
xcgfus = 12716.724 * 0.001
"Fuselage1 nose"
ycgfus1 = 0 * 0.001
xcgfus1 = 2.02
"Fuselage2 nose"
ycgfus2 = 0 * 0.001
xcgfus2 = (4.04 +(3.602/2))
"Tank cg locations"
ycgtank = 0 * 0.001
xcgtank = 18254 * 0.001
"Wing cg locations"
ycgwing = 0 * 0.001
xcgwing = 11010.341 * 0.001
"Horizontal tail cg locations"
ycght = 0 * 0.001
xcght = 22622.812 * 0.001
"Vertical tail cg locations"
ycgvt = 0 * 0.001
xcgvt = 21994.226 * 0.001
"Nacelle in. cg locations"
ycgnin = 4054.115 * 0.001
xcgnin = 8402.8547 * 0.001
"Nacelle out. cg locations"
ycgnout = 7100 * 0.001
xcgnout = 8927.0576 * 0.001
"Fuel Cell cg locations"
ycgfc = 0 * 0.001
xcgfc = 13000.91 * 0.001
"Main Land Gear cg locations"
ycgmlg = 1796.9 * 0.001
xcgmlg = 12000 * 0.001
"Nose Land Gear cg locations"
ycgnlg = 0 * 0.001
xcgnlg = 2470 * 0.001     #3511.67 * 0.001
"wing 1"
ycgwing1 = 0 * 0.001
xcgwing1 = 11010.341 * 0.001
"wing 2"
ycgwing2 = 2.096222
xcgwing2 = 11010.341 * 0.001
"wing 3"
ycgwing3 = 4.140828
xcgwing3 = 11010.341 * 0.001
"wing 4"
ycgwing4 = 6.083473
xcgwing4 = 11010.341 * 0.001
"wing 5"
ycgwing5 = 7.876322
xcgwing5 = 11010.341 * 0.001
"wing 6"
ycgwing6 = 9.475231
xcgwing6 = 11010.341 * 0.001
"wing 7"
ycgwing7 = 10.84083
xcgwing7 = 11010.341 * 0.001
"wing 8"
ycgwing8 = 11.93949
xcgwing8 = 11010.341 * 0.001
"wing 9"
ycgwing9 = 12.74416
xcgwing9 = 11010.341 * 0.001
"wing 10"
ycgwing10 = 13.23502
xcgwing10 = 11010.341 * 0.001

"Mass wing in parts"
Wwing1 = 2260.976941071543/9.81
Wwing2 = 2111.3460242596434/9.81
Wwing3 = 1909.7267479765471/9.81
Wwing4 = 1661.082981797737/9.81
Wwing5 = 1371.539912630189/9.81
Wwing6 = 1048.2254068128698/9.81
Wwing7 = 699.0963501006555/9.81
Wwing8 = 332.75508672282047/9.81
Wwing9 = 41.78655881294344/9.81
Wwing10 = 33.140695956115/9.81

SumW = Wwing1+Wwing2+Wwing3+Wwing4+Wwing5+Wwing6+Wwing7+Wwing8+Wwing9+Wwing10

lf = 23.876

"W are actually masses in [kg]"
W_fus = 0.174*m_oem
W_hor_tail = 0.048*m_oem
W_ver_tail =0.052*m_oem
Engine_weight = 0.04*m_oem
Fuel_Cell_Weight = 0.064*m_oem
LH2_system_tank = 0.075*m_oem
W_land_main = 0.041*m_oem
W_land_nose = 0.008*m_oem
W_controls = 0.09*m_oem
W_furnish = 0.083*m_oem
W_avionics = 0.04*m_oem
# mtot = W_fus+2*SumW+W_hor_tail+W_ver_tail+Engine_weight+Fuel_Cell_Weight+LH2_system_tank+W_land_main+W_land_nose
mtot =  W_avionics +W_fus+2*SumW+W_hor_tail+W_ver_tail+Engine_weight+Fuel_Cell_Weight+LH2_system_tank+W_land_main+W_land_nose
in_nacelle_percent  = 0.4
out_nacelle_percent = 0.1
W_fus1 = (4.04/lf)*W_fus
W_fus2 = (10.806/lf)*W_fus
W_fus3 = (9.03/lf)*W_fus
def Izz(mass,xcg,ycg):
    Izz = mass*((xcg-Xcg)**2+(ycg-Ycg))
    return Izz

xcg_control = xcgfus
ycg_control= 0
xcg_furnish = xcgfus
ycg_furnish= 0
xcg_avionics= xcgfus2
ycg_avionics= 0

Izz_fuse = Izz(W_fus,xcgfus, ycgfus)
Izz_tank = Izz(LH2_system_tank, xcgtank, ycgtank)
Izz_wing = Izz(SumW,xcgwing, ycgwing)
Izz_htail = Izz(W_hor_tail, xcght, ycght)
Izz_vtail = Izz( W_ver_tail, xcgvt, ycgvt)
Izz_nacelle_in = 2 * Izz((Engine_weight * in_nacelle_percent), xcgnin, ycgnin)
Izz_nacelle_out = 2 * Izz((Engine_weight * out_nacelle_percent), xcgnout, ycgnout)
Izz_fuel_cell = Izz((Fuel_Cell_Weight), xcgfc, ycgfc)
Izz_main_LG = 2 * Izz((W_land_main*0.5), xcgmlg, ycgmlg)
Izz_nose_LG = Izz((W_land_nose), xcgnlg, ycgnlg)
Izz_control = Izz(W_controls,xcg_control,ycg_control)
Izz_furnish = Izz(W_furnish,xcg_furnish,ycg_furnish)
Izz_Avionics=Izz(W_avionics,xcg_avionics,ycg_avionics)

total_Izz = Izz_fuse + Izz_tank + Izz_wing + Izz_htail + Izz_vtail + Izz_nacelle_in + Izz_nacelle_out + Izz_fuel_cell + Izz_main_LG +Izz_nose_LG
total_Izz = Izz_fuse + Izz_tank + Izz_wing + Izz_htail + Izz_vtail + Izz_nacelle_in + Izz_nacelle_out + Izz_fuel_cell + Izz_main_LG +Izz_nose_LG

kz_total =np.sqrt(total_Izz/mtot)

KZ = kz_total/(bw)
KZ2 = KZ**2

print("KZ^2 =", KZ2)
print("Izz =", total_Izz)

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

# print("KX2 = ",0.019 , "TBD")                          # squared Non-dimensional radius of gyration about the X-axis   [-]
# print("KY2 = ",1.25 * 1.114, "TBD")                    # squared Non-dimensional radius of gyration about the Y-axis   [-]
# print("KZ2 = ",0.042 , "TBD")                          # squared Non-dimensional radius of gyration about the Z-axis   [-]
# print("KXZ = ",0.002  , "TBD")                         # Non-dimensional product of inertia                            [-]

print("CXde =", Delta_e_Derivatives()[0])
print("CZde =", Delta_e_Derivatives()[1])
print("Cmde =", Delta_e_Derivatives()[2])