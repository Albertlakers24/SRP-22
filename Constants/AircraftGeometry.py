''' AIRCRAFT GEOMETRY OUTPUTS --> GABRIEL  ---> Last Update - 27/2/22 '''
import numpy as np

# Suggested headers - pls change if needed -- segregate as you'd like

'''TO BE UPDATED: High Lift Devices (for Alberto)''' #todo to be updated!
twist = 4/180*np.pi             # Wing twist angle                          [rad]
dihedralw = 0                   # Wing dihedral angle                       [rad]

DCL_Alpha0_TO =-6.00019144402   # Change Alpha_0 for take off               [deg]
DCL_Alpha0_LD =-9.0002871660244 # Change Alpha_0 for landing                [deg]
delta_CL_TO = 0.49626601583949  # delta CL for take-off                     [-]
delta_CL_LD = 0.853914517995425 # delta CL for landing                      [-]
CL_Alpha_TO = 6.204527397612    # change of CL_Alpha take-off               [-]
CL_Alpha_LD = 6.365986907313    # change of CL_Alpha landing                [-]


''' Wing planform'''
taperw = 0.45                   # Wing taper ratio                                  [-]
S_w = 65.28286905               # Wing surface area                                 [m^2]
Aw   = 12                       # Wing aspect ratio                                 [-]
bw   = 27.9891841360194        # Wing span                                         [m]
c_tw = 1.4477164208285898          # Wing tip chord                                    [m]
c_rw = 3.2171476018413103        # Wing root chord                                   [m]
t_c_ratio_w = 0.21      # Wing thickness to chord ratio                     [-]
c_mac_w =   2.4442926032380528     # Length of wing mean aerodynamic chord             [m]
y_mac   =   6.112580443498489       # Spanwise position of wing mean aerodynamic chord  [m]
Sweep_quarterchordw = 0 # Sweep at quarter chord of the wing                [rad]   todo: update the drawing with these values later
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
l_nc = 3.7              # length of nosecone [m]
n_row = 12              # Number of rows of Seat
skin_thickness = 0.105  # Fuselage skin thickness [m]

w_door_front = 0.61
w_door_back = 0.508
l_lav = 0.914           # lavatory length   [m]
w_lav = 0.914           # lavatory width    [m]
l_galley = 0.762        # galley length     [m]
w_galley = 0.914        # galley width      [m]
l_seat = 0.76           # seat pitch        [m]
l_tank = 0              #
l_pax = l_seat * n_row

'''High Lift Devices and Aileron'''
SwfS_flap = 0.600451            # Wing area by flaps area (Swf/S)           [-]
S_flap = 8.48126846893824       # Flaps area                                [m^2]
c_accent_c_flap = 1.216364      # c'/c for the flaps                        [-]
b_flap = 7.38 - 1.5             # Span flaps                                [m]
cf_over_cprime_flap = 0.2877429 # Chord flap/c'                             [-]
cf_cmac = 0.35                  # Chord flap/ MAC                           [-]
delta_a_left = 20.0375           # Left deflection aileron                   [deg]
delta_a_right = 26.17143        # Right deflection aileron                  [deg]
cf_over_c_aileron = 0.30        # chord flap/c for the aileron              [-]
y_outboard_aileron = 12.587     # y outboard location of the aileron        [m]
y_inboard_aileron = 7.487       # y inboard location of the aileron         [m]
y_outboard_flap = 7.38          # y outboard location of the flap           [m]
y_inboard_flap = 1.5            # y inboard location of the flap            [m]
delta_flap =40/180*np.pi        # Flap deflection                           [rad]