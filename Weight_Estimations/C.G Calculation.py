import numpy as np
from Constants.AircraftGeometry import l_f,d_f_outer, c_mac_w, l_nc, w_door_front, l_pax, w_lav, l_tank
# Some parameters
w_f = d_f_outer                # Width of fuselage
h_f = w_f                    # Height of fuselage, assumed to be a circular fuselage
x_cg_LEMAC = 0.35 * c_mac_w  # (predetermined value)
#distance of fuselage group components w.r.t. datum
x_fuselage = 0.4 * l_f
x_empennage = 0.9 * l_f
x_sys = 0.4 * l_f
x_wing = 0.4 * c_mac_w
x_prop_wing = -0.2 * c_mac_w
x_prop_fuse = 0.85 * l_f
x_hydrogen_tank = l_nc + w_door_front + l_pax + w_lav + l_tank/2
Mf_fuselage = 0.12  # Mass fraction of components w.r.t MTOM
Mf_prop = 0.07
Mf_empennage = 0.07
Mf_sys = 0.18
Mf_wing = 0.11
Mf_tank = 0.09    # Mass fraction of the tank to be updated
print("======================================")
print("For hydrogen fuel cell architecture: ")
print("hydrogen tank placement", x_hydrogen_tank/l_f)
M_fg_frac = Mf_fuselage + Mf_empennage + Mf_sys + Mf_tank  # sum of mass of the fuselage group
M_wg_frac = Mf_wing + Mf_prop  # sum of mass of the wing group
x_fg = (Mf_fuselage * x_fuselage + Mf_empennage * x_empennage + Mf_sys * x_sys + Mf_tank * x_hydrogen_tank) / M_fg_frac  # c.g. of the fuselage group
x_wg_LEMAC = (Mf_wing * x_wing + Mf_prop * x_prop_wing) / M_wg_frac
x_LEMAC = x_fg - x_cg_LEMAC + M_wg_frac / M_fg_frac * (x_wg_LEMAC - x_cg_LEMAC) # The distance from nose to LEMAC.
x_oewcg = x_LEMAC + x_cg_LEMAC        # The distance from nose to cg.

print("The distance from zero point to LEMAC is ",x_LEMAC, "[m]")
print("The distance from zero point to cg is ", x_oewcg, "[m]")