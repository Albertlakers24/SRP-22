''' EMPENNAGE & LANDING GEAR OUTPUTS --> INES  ---> Last Update - 13/3/22 '''
import numpy as np
# Combined the two because I thought they weren't interrelated and I wanted to make fewer files

tc_tail = 0.09                     # Average t/c ratio of wing crosssection
xc_m_tail = 0.29                        # (x/c)_m Position of max thickness - 0.3 for low speed, 0.5 for high speed

'''To be updated :'''
lambdahalf_h = 0                # Sweep at half chord of HT     [rad]
cmac_h =  1.614 #1              # MAC for HT                    [m]
bh = 6                          # Span of the HT                [m]
c_rh = 1.98 #2                  # HT root chord                 [m]
c_th =  1.29 #1                 # HT tip chord                  [m]
cmac_v =  3.204                 # MAC for VT                    [m]
c_rv = 3.88                     # VT  root chord                [m]
c_tv = 1.98                     # VT tip chord                  [m]

''' Vertical Tail '''
Sv = 13.96                      # Surface area VT               [m^2]
taperv = 0.5                    # Taper VT                      [-]
Vh_V = 1                        # Volume fraction for T-tail    [-]
#Vh_V = np.sqrt(0.85)                    # Volume fraction for conventional tail [-]

''' Horizontal tail '''
Sh = 11.376                     # Surface area HT               [m^2]
taperh = 0.75                   # Taper HT                      [-]
x_h = 23.3                      # x-location horizontal tail    [m]
A_h = 0                         # HT aspect ratio               [-]
lh = 10.6                       # Tail arm                      [m]


''' Landing Gear '''
