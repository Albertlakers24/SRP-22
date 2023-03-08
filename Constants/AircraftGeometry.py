''' AIRCRAFT GEOMETRY OUTPUTS --> GABRIEL  ---> Last Update - 27/2/22 '''

# Suggested headers - pls change if needed -- segregate as you'd like


''' Wing planform'''
taper = 0.45            # wing taper ratio                                  [-]
S_w = 59.8947           # Wing surface area                                 [m^2]
Aw   = 12               # Wing aspect ratio                                 [-]
bw   = 26.81071         # Wing span                                         [m]
c_t = 1.38661           # Wing tip chord                                    [m]
c_r = 3.08136           # Wing root chord                                   [m]
t_c_ratio = 0.21        # Wing thickness to chord ratio                     [-]
c_mac   =   2.34113     # length of wing mean aerodynamic chord             [m]
y_mac   =   5.855       # spanwise position of wing mean aerodynamic chord  [m]


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