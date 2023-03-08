import numpy as np
import scipy as sp
from Initial_Aircraft_Sizing.Wing_planform import taper
#Switch
switch = 1          # put 1 for circular fuselage, 2 for double bubble fuselage
prop_choice = 4 # 1 = H2 combustion, 2 = series/ parallel, 3 = series, 4 = H2 fuel cell

# Some constants
k_cabin = 1.08      #For single aisle

# Dimensions under the regulation
w_door_front = 0.61
w_door_back = 0.508
l_lav = 0.914
w_lav = 0.914
l_galley = 0.762
w_galley = 0.914
l_seat = 0.76      # seat pitch [m]

# Design choices
n_SA = 4         # Number of seats abreast

if n_SA == 3:
    n_PAX = 48                  # Number of passengers (46-50)
    D_inner = 2.48              # inner diameter of the fuselage cross section [m]
elif n_SA == 4:
    n_PAX = 48                  # Number of passengers (46-50)
    D_inner = 2.8               # inner diameter of the fuselage cross section [m]
elif n_SA == 5:
    n_PAX = 50                  # Number of passengers (46-50)
    D_inner = 3.5               # inner diameter of the fuselage cross section [m]
n_row = np.ceil(n_PAX/n_SA)     # number of rows
l_cp = 4                        # cockpit length [m]
double_bubble_height = 2.7

if switch == 1:
    height = D_inner
else:
    height = double_bubble_height

# Basic Length/Dimensions calculations
D_eff = np.sqrt(height*D_inner)             # effective (inner) diameter [m]
l_cabin = n_row * k_cabin                   # Length of the cabin [m]
skin_thickness = 0.084 + 0.045 * D_eff      # fuselage skin thickness [m]
D_outer = D_eff + skin_thickness            # outer fuselage diameter [m]




if prop_choice == 1:
    if n_SA == 3:
        l_tank = 6.4
    if n_SA == 4:
        l_tank = 5.1                            # Length of the fuel tank as fuselage section
    if n_SA == 5:
        l_tank = 3.26
elif prop_choice == 4:
    if n_SA == 3:
        l_tank = 2.59
    if n_SA == 4:
        l_tank = 2.1                            # Length of the fuel tank as fuselage section
    if n_SA == 5:
        l_tank = 1.3446

else:
    l_tank = 0

l_t = 1.6 * D_outer + l_tank                      # Length of tail [m]
l_f = l_cabin + l_t + l_cp +1.02                       # length of the fuselage [m]
fineness_ratio = l_f/D_outer                      # Fineness ratio
l_nc = 1.2 * D_outer                                # length of the nose cone [m]
l_tc = 3 * D_outer                                # length of the tail cone [m]
l_pax = l_seat * n_row                            # length of the passenger seating area [m]
l_constant_cross_section = l_f - l_nc - l_tc      # length of the cross section with constant cross section [m]
theta_tc = np.arctan(D_outer/l_tc)*180/np.pi
print(l_f)
# Skin Drag Computation



# Fuselage mass estimation (Raymer p475)
K_ws = 0.75*(1+2*taper)/(1+taper)
K_door = 1.12   # 1.0 if no cargo door, 1.06 if one side cargo door,
                # 1.12 if two side cargo door, 1.12 if aft clamshell door
                # 1.25 if two side cargo doors and aft clamshell door
K_lg = 1     # 1.12 if fuselage mounted, 1.0 if otherwise
MTOM = 30000
W_dg = MTOM * 9.80665 * 0.224809      # Flight design gross weight [lb] "The aircraft weight at which the structure will withstand the design load factor”
S_f_wet = l_constant_cross_section * D_outer + D_outer * 0.5 * l_nc + D_outer * 0.5 * l_tc                            # Fuselage wetted area [m^2]
S_f_wet_imp = (l_constant_cross_section*D_outer+D_outer*0.5*l_nc+D_outer*0.5*l_tc) * 10.7639            # Fuselage wetted area in imperial [ft^2]
L = (l_f - l_t) * 3.28084     # Fuselage structural length (excl. tailcap) [ft]
D = D_eff * 3.28084           # Structural Depth [ft]
N_Z = 1.5*2.5                 # 1.5 * limit load factor
W_fuselage = 0.3280 * K_door * K_lg * ((W_dg*N_Z)**0.5) * (L**0.25)* (S_f_wet_imp**0.302)* ((L/D)**0.1) * ((1+K_ws)**0.04)      # Weight of the fuselage [lb]
W_fuselage_SI = W_fuselage * 4.44822    # Weight of the fuselage [N]
Mass_fuselage = W_fuselage * 4.44822 / 9.80665    # mass` of the fuselage [kg]