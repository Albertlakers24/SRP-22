''' AIRCRAFT GEOMETRY OUTPUTS --> GABRIEL  ---> Last Update - 27/2/22 '''
import numpy as np

# Suggested headers - pls change if needed -- segregate as you'd like

'''TO BE UPDATED: High Lift Devices (for Alberto)''' #todo to be updated!
SwfS_flap = 1                   # Wing area by flaps area (Swf/S)           [-]
c_accent_c_flap = 1             # c'/c for the flaps                        [-]
b_flap = 1                      # Span flaps                                [m]
cf_over_cprime_flap = 1         # Chord flap/c'                             [-]
delta_a_left = 1                # Left deflection aileron                   [rad]
delta_a_righ = 1                # Right deflection aileron                  [rad]
cf_over_c_aileron = 0.38        # chord flap/c for the aileron              [-]
y_outboard_aileron = 7          # y outboard location of the aileron        [m]
y_inboard_aileron = 1           # y inboard location of the aileron         [m]
y_outboard_flap = 1             # y outboard location of the flap           [m]
y_inboard_flap = 1              # y inboard location of the flap            [m]
delta_flap =40/180*np.pi        # Flap deflection                           [rad]
twist = 0                       # Wing twist angle                          [rad]
dihedralw = 0                   # Wing dihedral angle                       [rad]

''' Wing planform'''
taperw = 0.45           # wing taper ratio                                  [-]
S_w = 63.6872413973     # Wing surface area                                 [m^2]
Aw   = 12               # Wing aspect ratio                                 [-]
bw   = 27.645015        # Wing span                                         [m]
c_tw = 1.38661          # Wing tip chord                                    [m]
c_rw = 3.08136          # Wing root chord                                   [m]
t_c_ratio_w = 0.21      # Wing thickness to chord ratio                     [-]
c_mac_w =   2.34113     # Length of wing mean aerodynamic chord             [m]
y_mac   =   5.855       # Spanwise position of wing mean aerodynamic chord  [m]
Sweep_quarterchordw = 0 # Sweep at quarter chord of the wing                [rad]   todo: check
Sweep_halfchordw = -1.81*np.pi/180 # Sweep at half chord of the wing        [rad]
SweepLE = 1.81*np.pi/180 # Sweep at LE of the wing                          [rad]

''' Fuselage geometry'''
# 4SA configuration
l_f     = 23.964        # Length of the fuselage    [m]
l_tc    = 9.01909       # Length of the tailcone    [m]
d_f_inner = 2.8         # fuselage inner diameter
d_f_outer = 3.01        # fuselage outer diameter
slenderness = 7.96146   # fuselage slenderness ratio
l_t     = 8.46208       # fuselage tail length [m]
                        # length behind the aft pressure bulkhead
l_cp            = 3.71091       # Length of the cockpit [m]
l_cylindrical   = 11.24491      # Length of the cylindrical section of the fuselage [m]
l_cabin = 11.79101      # length of cabin [m]
n_row = 12              # Number of rows of Seat
skin_thickness = 0.105  # Fuselage skin thickness [m]

w_door_front = 0.61
w_door_back = 0.508
l_lav = 0.914           # lavatory length   [m]
w_lav = 0.914           # lavatory width    [m]
l_galley = 0.762        # galley length     [m]
w_galley = 0.914        # galley width      [m]
l_seat = 0.76           # seat pitch        [m]
