import numpy as np
#from Initial_Aircraft_Sizing.Wing_planform import Sw, c_mac,specific_gas_constant, gamma, M_cross,b
#from Initial_Aircraft_Sizing.Fuselage import l_f


##Design 1; Fuel Cell
Design = 1

# Inputs: Constants from other files
D_inner = 3.0                       # Inner diameter        [m]
thickness = 0.21                    # Fuselage thickness    [m]
D_outer = D_inner+ thickness*2      # Outer diameter        [m]
bw = 26.8                           # Wing span             [m]
Sw = 59.9                           # Wing surface area     [m^2]
l_f = 23.9                          # Fuselage lenght       [m]
c_mac = 2.3                         # MAC wing              [m]
xcg_aft = 12.7                      # Aft cg location       [m]
l_ultimate = 19.755                 # Ultimate location     [m]     TECHNICAL DRAWING INPUT

if Design ==1:
    # Inputs: Variables
    x_h = 23.3                      # Location horizontal tail  [m]
    x_v = 21.9                      # Location vertical tail    [m]
    Vh = 0.86                       # Horizontal tail volume    [m^3]   -> 0.68113 (propeller)
    Vv = 0.08                       # Vertical tail volume      [m^3]   -> 0.08 (expected for propeller)
    Av = 1.3                        # VT aspect ratio           [-]     -> 1-2
    Ah = 4                          # HT aspect ratio           [-]     -> 3-5
    taperv = 0.5                    # VT taper ratio            [-]     -> 0.3-0.7
    taperh = 0.75                   # HT taper ratio            [-]     -> 0.3-1
    Kc = 1.4                        # ?

    # Calculations
    l_v = x_v - xcg_aft             # VT location               [m]
    l_h = x_h - xcg_aft             # HT location               [m]
    l_opt = Kc * np.sqrt((4 * c_mac * Sw * Vh) / (np.pi * D_outer))

    print("-----------Location check---------")
    print("l_v =", l_v)
    print("l_opt=", l_opt)

# HORIZONTAL TAIL
# Calculations changes for Conventional =1, Conventional plus dorsal =2
switch =1
Type =1
#abs((Sh-Sh_stability)/Sh)*100 < 10:

if Type ==1:
    Sh = Vh * (Sw * c_mac) / l_h               # HT surface area           [m^2]
    if switch ==1:

        # Determine Horizontal Tail Geometry
        bh = np.sqrt(Sh*Ah)                 # HT wing span              [m]
        c_rh = (2*Sh)/((1+taperh)*bh)       # HT root chord             [m]
        c_th = c_rh * taperh                # HT tip chord              [m]
        t_c_ratioh = 0.18                   # t\c -> ASSUMPTION min(0.18, ((M_cross - M_dd) - 0.115 * (C_Lhat ** 1.5)))        # no sweep in horizontal tail! thickness to chord ratio
        c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))     # HT length of MAC          [m]
        y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh             # Spanwise location of MAC  [m]

        # PRINT STATEMENT HORIZONTAL TAIL INITIAL GEOMETRY
        print("-----------Horizontal Tail Sizing------------")
        print("Sh =", Sh, "m^2")
        print("bh =", bh, "m")
        print("cr_h =", c_rh, "m")
        print("ct_h =",c_th, "m")
        print("c_mach =", c_mach_h)
        print("y_machh =", y_mach_h)
        print("ratio wing areas =", Sh/Sw)
        print("tail_arm_h=", x_h - xcg_aft)

        Sv = Vv * (Sw * bw) / l_v           # VT surface area           [m^2]
        bv = np.sqrt(Av*Sv)                 # VT wing span              [m]
        c_mac_v = Sv/bv                     # VT MAC                    [m]
        c_rv = 2*Sv/((1+taperv)*bv)         # VT root chord             [m]
        c_tv = taperv*c_rv                  # VT tip chord              [m]
        y_mach_v = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv # [m]
        #c_mac_v =(2/3)*c_rv*((1+taperv+taperv**2)/(taperv+1))

        # PRINT STATEMENT VERTICAL TAIL INITIAL GEOMETRY
        print("-----------Vertical Tail Sizing------------")
        print("Sv =", Sv, "m^2")
        print("bv =", bv, "m")
        print("cr_v =", c_rv, "m")
        print("ct_v =",c_tv, "m")
        print("c_macv =", c_mac_v)
        print("y_machv =", y_mach_v)
        print("tail_arm_v=", x_v-xcg_aft)

        # print("Sh ratio", Sw/Sh)
        # print("Sv ratio", Sw/Sv)

        # print("----------------DRAWING CALCULATIONS-------------")
        # print("Location Horizontal Tail =", xcg_aft + l_opth)
        # print("Distance from TE to tail =", l_f - (xcg_aft + l_opth))
    #     print("l_h", l_h)
    #     print("l_v", l_v)
    #     print("part outside", (x_v+(3/4)*c_rv)-l_f)
    #     print("x_h=", x_h)
    #     print("x_v =", x_v)
    #     print("nose to LE vertical tail=", x_v-(1/4*c_rv))
    #     print("nose to LE horizontal tail=", x_h-(1/4*c_rh))
    #     if (x_v-(1/4)*c_rv)-l_ultimate < 0:
    #         print("we are interfering :( by", (x_v - (1 / 4) * c_rv) - l_ultimate)
    #     else:
    #         print("not interfering by:",(x_v - (1 / 4) * c_rv) - l_ultimate)
    # else:
    #     print("Sh =", Sh, "m^2")
    #     print("Rerun Code for lh")



#print("distance outside fuselage=", l_f-(xcg+l_opt+c_rh))
#Check with graph
#Sh_stability = 11.7                 # Determine from the stability graph
#CHOOSE
#lh = 24                            # m from CG calculations I think
#lv = 24                            # m from CG calculations
#sweep_df_LE = 74                    # degrees TBD
#sweeph_quarter = 0                  # degrees (0 for propeller)
#sweepv_LE = 25                      # degrees (0-50) propeller
#sweepv_half = 25                    # degrees -> depending on the t/c
#fraction_df = 0.19                  # -






"""
if Type ==2:
    Sh = Vh * (Sw * cw) / lh  # m^2 Horizontal tail surface area
    if abs((Sh - Sh_stability) / Sh) * 100 < 10:
        # Determine span, root chord, tipchord MAC
        bh = np.sqrt(Sh * Ah)  # m horizontal tail span
        c_rh = (2 * Sh) / ((1 + taperh) * bh)  # m root chord
        c_th = c_rh * taperh  # m Root chord
        t_c_ratioh = min(0.18, ((M_cross - M_dd) - 0.115 * (
                    C_Lhat ** 1.5)))  # no sweep in horizontal tail! thickness to chord ratio
        c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))  # length of MAC
        y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh  # Spanwise location of MAC
        #Size Dorsal
        S_df  ##Calculate!!!
        h_df = np.sqrt((2 * S_df) / (np.tan(sweep_df_LE) - np.tan(sweepv_LE)))
        # Position the surface
        # Print Statements
        print("Sh =", Sh, "m^2")
if Type ==3:
    print("H")
#VERTICAL TAIL
if Type ==1:
    #Choose longitudinal location of vertical tail aerodynamic center
    #determine moment arm to the most aft CG position - lv
    #Compute required tail area:
    Sv = Vv*(Sw*bw)/lv                      # m^2 Horizontal tail surface area
    #Determine sweep, AR, taper
    #Determine span, cr, ct, MAC
    bv = np.sqrt(Sv * Av)  # m horizontal tail span
    c_rv = (2 * Sv) / ((1 + taperv) * bv)  # m root chord
    c_tv = c_rv * taperv  # m Root chord
    t_c_ratiov = min(0.18, (((np.cos(sweepv_half))**3)*(M_cross - M_dd) - 0.115 * (C_Lhat ** 1.5))/(np.cos(sweepv_half))**2)  # no sweep in vertical tail! thickness to chord ratio
    c_machv = (2 / 3) * c_rv * ((1 + taperv + taperv ** 2) / (1 + taperv))  # length of MAC
    y_machv = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv
    #position the surface
if Type ==2:
    print("goodbye")
if Type ==3:
    print("H")"""