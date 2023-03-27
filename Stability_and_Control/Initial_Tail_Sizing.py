import numpy as np
from Constants.AircraftGeometry import d_f_outer, bw, S_w, c_mac_w
from Constants.Masses_Locations import xcg_aft_potato

print("FILE: Initial Tail Sizing")

Initial = True

# Inputs: Esimated Variables
x_h = 23.3                      # Location horizontal tail  [m]
x_v = 21.9                      # Location vertical tail    [m]
Vh = 0.86                       # Horizontal tail volume    [m^3]   -> 0.68113 (propeller)
Vv = 0.08                       # Vertical tail volume      [m^3]   -> 0.08    (expected for propeller)
Av = 1.3                        # VT aspect ratio           [-]     -> 1-2
Ah = 4                          # HT aspect ratio           [-]     -> 3-5
taperv = 0.5                    # VT taper ratio            [-]     -> 0.3-0.7
taperh = 0.75                   # HT taper ratio            [-]     -> 0.3-1
Kc = 1.4
Sh_newvalue = 11.376            # HT surface area           [m^2]
Sv_newvalue = 13.96             # VT surface area           [m^2]

# Intermediate Calculations
l_v = x_v - xcg_aft_potato      # VT location               [m]
l_h = x_h - xcg_aft_potato      # HT location               [m]
l_opt = Kc * np.sqrt((4 * c_mac_w * S_w * Vh) / (np.pi * d_f_outer))

print("-----------Location check---------")
print("l_v =", l_v)
print("l_opt=", l_opt)

if Initial == True:
    Sh = Vh * (S_w * c_mac_w) / l_h
    Sv = Vv * (S_w * bw) / l_v
else:
    Sh = Sh_newvalue
    Sv = Sv_newvalue

# HORIZONTAL TAIL
def Geometry_HT():
    bh = np.sqrt(Sh * Ah)
    c_rh = (2 * Sh) / ((1 + taperh) * bh)
    c_th = c_rh * taperh
    t_c_ratioh = 0.18
    c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))
    y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh
    return bh, c_rh,c_th,c_mach_h,y_mach_h

def Geometry_VT():
    bv = np.sqrt(Av * Sv)
    c_mac_v = Sv / bv
    c_rv = 2 * Sv / ((1 + taperv) * bv)
    c_tv = taperv * c_rv
    c_mac_v = (2 / 3) * c_rv * ((1 + taperv + taperv ** 2) / (taperv + 1))
    y_mach_v = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv
    return bv,c_rv,c_tv,c_mac_v,y_mach_v

print("-----------Horizontal Tail Sizing------------")
print("Sh =", Sh, "m^2")
print("bh =", Geometry_HT()[0], "m")
print("cr_h =", Geometry_HT()[1], "m")
print("ct_h =", Geometry_HT()[2], "m")
print("c_mach =", Geometry_HT()[3],"m")
print("y_machh =", Geometry_HT()[4],"m")
print("ratio wing areas =", Sh / S_w)
print("tail_arm_h=", x_h - xcg_aft_potato, "m")

print("-----------Horizontal Tail Sizing------------")
print("Sh =", Sv, "m^2")
print("bh =", Geometry_VT()[0], "m")
print("cr_h =", Geometry_VT()[1], "m")
print("ct_h =", Geometry_VT()[2], "m")
print("c_mach =", Geometry_VT()[3],"m")
print("y_machh =", Geometry_VT()[4],"m")
print("ratio wing areas =", Sv / S_w)
print("tail_arm_h=", x_v - xcg_aft_potato, "m")