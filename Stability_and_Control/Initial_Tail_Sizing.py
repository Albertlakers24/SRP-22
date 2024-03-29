import numpy as np
from Constants.AircraftGeometry import d_f_outer, bw, S_w, c_mac_w
from Constants.Masses_Locations import xcg_aft_potato

print("FILE: Initial Tail Sizing")

Initial = True

# Inputs: Esimated Variables
x_h = 21.45                     # HT: nose to ac of HT      [m]
x_v = 20.9                      # VT: nose to ac of VT      [m]
Vh = 0.86                       # Horizontal tail volume    [m^3]   -> 0.68113 (propeller)
Vv = 0.08                       # Vertical tail volume      [m^3]   -> 0.08    (expected for propeller)
Av = 1.3                        # VT aspect ratio           [-]     -> 1-2
Ah = 4                          # HT aspect ratio           [-]     -> 3-5
taperv = 0.5                    # VT taper ratio            [-]     -> 0.3-0.7
taperh = 0.75                   # HT taper ratio            [-]     -> 0.3-1
Kc = 1.4
Sh_newvalue = 11.376            # HT surface area           [m^2]
Sv_newvalue = 13.96             # VT surface area           [m^2]
Sweepquarterchord_h = 0         # HT c/4 sweep              [rad]
Sweep_LE_vt = 30*np.pi/180       # VT LE sweep               [rad]

# Intermediate Calculations
l_v = x_v - xcg_aft_potato      # VT location               [m]
l_h = x_h - xcg_aft_potato      # HT location               [m]
l_opt = Kc * np.sqrt((4 * c_mac_w * S_w * Vh) / (np.pi * d_f_outer))

print("-----------Location check---------")
print("x_v =", x_v)
print("x_h =", x_h)
print("l_v =",l_v)
print("l_h = ", l_h)

if Initial == True:
    Sh = Vh * (S_w * c_mac_w) / l_h
    Sv = Vv * (S_w * bw) / l_v
else:
    Sh = Sh_newvalue
    Sv = Sv_newvalue

# HORIZONTAL TAIL
def Geometry_HT(Sh, Ah, taperh):
    bh = np.sqrt(Sh * Ah)
    c_rh = (2 * Sh) / ((1 + taperh) * bh)
    c_th = c_rh * taperh
    t_c_ratioh = 0.18
    c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))
    y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh
    d_LE_h = x_h - c_mach_h/4                                               # Distance nose a/c to LE HT
    return bh, c_rh,c_th,c_mach_h,y_mach_h,d_LE_h

def Geometry_VT():
    bv = np.sqrt(Av * Sv)
    c_mac_v = Sv / bv
    c_rv = 2 * Sv / ((1 + taperv) * bv)
    c_tv = taperv * c_rv
    c_mac_v = (2 / 3) * c_rv * ((1 + taperv + taperv ** 2) / (taperv + 1))
    y_mach_v = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv
    d_LE_v = x_v - c_mac_v/4 - y_mach_v*np.tan(Sweep_LE_vt)                     # Distance nose a/c to LE VT
    return bv,c_rv,c_tv,c_mac_v,y_mach_v,d_LE_v

print("-----------Horizontal Tail Sizing------------")
print("d_LE_h =", Geometry_HT(Sh, Ah, taperh)[5], "m")
print("Sh =", Sh, "m^2")
print("Ah =", Ah)
print("bh =", Geometry_HT(Sh, Ah, taperh)[0], "m")
print("cr_h =", Geometry_HT(Sh, Ah, taperh)[1], "m")
print("ct_h =", Geometry_HT(Sh, Ah, taperh)[2], "m")
print("c_mach =", Geometry_HT(Sh, Ah, taperh)[3],"m")
print("y_machh =", Geometry_HT(Sh, Ah, taperh)[4],"m")
print("ratio wing areas =", Sh / S_w)
print("tail_arm_h=", x_h - xcg_aft_potato, "m")

print("-----------Vertical Tail Sizing------------")
print("d_LE_v =", Geometry_VT()[5], "m")
print("Sv =", Sv, "m^2")
print("Av =", Av)
print("bv =", Geometry_VT()[0], "m")
print("cr_v =", Geometry_VT()[1], "m")
print("ct_v =", Geometry_VT()[2], "m")
print("c_macv =", Geometry_VT()[3],"m")
print("y_machv =", Geometry_VT()[4],"m")
print("ratio wing areas_v =", Sv / S_w)
print("tail_arm_v=", x_v - xcg_aft_potato, "m")