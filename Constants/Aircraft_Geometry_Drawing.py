''' Aircraft Geometry Drawing Outputs --> Gabriel  ---> Last Update - N/A '''

"""
Geometry calculated using the drawing:
* Distances
* Area hall?
"""

"""To be updated/added:"""
center_S_BS = 11            # Center location Side area aircraft    [m]



ln1 = 2.19543   # [m] Longitudinal distance between the big prop to the wing c/4
ln2 = 1.88328   # [m] Longitudinal distance between the small prop to the wing c/4
bn1 = 0.77      # [m] Big nacelle width (in y direction)
bn2 = 0.51      # [m] Small nacelle width (in y direction)
v_t_w = 3.83079 # [m] Verticl distance (z-direction) from the ac of the HT to the ac of wing (always +ve)
y_T = 4.05412   # [m] Lateral distance (y-direction) between the fuselage center line to the big engine.
Zw = -1.82253   # [m] Vertical distance (z-direction) between the fuselage
                # center line to the LE of the wing, assuming no dihedral.
l_fn = 10.24    # [m] Length nose to LE wing
zv = 2.422      # [m] Distance: cg aircraft to the ac of the VT (positive when going up)
zw = -1.84      # [m] Center line fuselage to root of the wing in the z direction -> negative for high wing

# center line to the LE of the wing, assuming no dihedral.
"""Area"""
S_BS = 74.8936  # Body side area [m**2]


''' updated 20-03 - fuselage '''
Amax_fuse = 7.1158                      # Max frontal area of fuselage [m^2]
A_base_fuse = 0.7 #1.138525             # Base area of fuselage (from the back of the fuselage) [m^2]
l_nacelle = 1.1975                           # Length of nacelle [m] (can still alter if we decided to install the MLG in the nacelles)
Amax_nacelle = 0.74375                     # Max frontal area of nacelle [m^2]
upsweep_ang =  3.03116                    # Upsweep angle of fuselage [deg] (centerline of fuselage to center of base - adsee II, lec 2, pg 51)


''' will be updated once the design is finalized - undercarriage '''
depth_lg = (0.838*4) + (0.343*2)                           # Depth of landing gear (strut + wheel)
width_lg =(0.247 *4) + (0.108*2)                           # Width of landing gear (strut + wheel)
S_lg_front = (0.834*0.247)*4 + (0.343*0.108*2)             # Frontal area of landing gear (strut area ) --- p.52 of ADSEE2, Lecture 2


# '''Distances'''
# ln1 = 1.983             # Big engine: front nacelle to c/4 of the wing      [m]
# ln2 = 1.34              # Small engine: front nacelle to c/4 of the wing    [m]
# bn1 = 0.766             # Big engine: width of nacelle                      [m]
# bn2 = 0.596             # Small engine: width of nacelle                    [m]
# v_t_w = 4.166           # Vertical distance HT to wing                      [m] todo: make more exact
# zv = 3.3026             # Center line fuselage to aerodynamic center VT     [m]
# y_T = 4                 # Center line fuselage to engine                    [m]
# Zw = -1.84              # Center line fuselage to wing                      [m] -> negative for high wing
# center_S_BS = 11        # x center location of the side area                [m]
# #l_nc    = 3.7           # Length of the nose cone                           [m]
#
# '''Area'''
# S_BS = 63.89            # Bode side area                                    [m^2]
