''' EMPENNAGE & LANDING GEAR OUTPUTS --> INES  ---> Last Update - 13/3/22 '''
import numpy as np

'''To be updated:'''
Av_eff = 3.5                    # VT: effective Aspect Ratio    [-]         FILE: Sizing_Vertical_Tail
xc_m_tail = 0.29                # (x/c)_m Position of max thickness - 0.3 for low speed, 0.5 for high speed
t_over_c_HT = 0.12              # HT : t/c                      [-]


''''First Fidelity: TAILARM'''
Vh_V = 1                        # Volume fraction for T-tail    [-]
Sh = 15.506                     # HT: Surface area              [m^2]
taperh = 0.75                   # HT: Taper                     [-]
x_h = 21.45                     # HT: xlocation from nose to ac [m]
A_h = 4                         # HT: aspect ratio              [-]
lh = 8.85                       # HT: Tail arm                  [m]
Sv = 19.33                      # VT: Surface area              [m^2]
taperv = 0.5                    # VT: Taper                     [-]
Av = 1.3                        # VT: Aspect ratio              [-]
tc_tail_vt = 0.12               # VT: t/c ratio                 [-]
Sweep_LE_HT = 3.94*np.pi/180    # HT: LE sweep angle            [rad]
Sweep_LE_VT = 30*np.pi/180      # VT: LE sweep angle            [rad]
Sweep_quarter_chord_HT = -2.73*np.pi/180 # HT: c/4 sweep angle  [rad]
Sweep_halfchord_VT = 16.1*np.pi/180     # VT: c/2 sweep angle   [rad]
lambdahalf_h = 0                # HT: c/2 sweep angle           [rad]
cmac_h =  1.982                 # HT: MAC                       [m]
bh = 7.8756                     # HT: Span                      [m]
c_rh = 2.2502                   # HT: root chord                [m]
c_th =  1.6876                  # HT: tip chord                 [m]
cmac_v =  3.817                 # VT: MAC                       [m]
c_rv = 4.908                    # VT: root chord                [m]
c_tv = 2.45                     # VT: tip chord                 [m]
lv = 8.30                       # VT: tail arm                  [m]
bv = 4.78                       # VT: Span                      [m]

''' Landing Gear '''

'''Old Values'''
# Sweep_quarter_chord_H = 0       # HT c/4 sweep angle            [rad]
# Sweep_halfchord_VT = 0          # VT c/2 sweep angle            [rad]
# lambdahalf_h = 0                # Sweep at half chord of HT     [rad]
# cmac_h =  1.614 #1              # MAC for HT                    [m]
# bh = 6                          # Span of the HT                [m]
# c_rh = 1.98 #2                  # HT root chord                 [m]
# c_th =  1.29 #1                 # HT tip chord                  [m]
# cmac_v =  3.204                 # MAC for VT                    [m]
# c_rv = 3.88                     # VT  root chord                [m]
# c_tv = 1.98                     # VT tip chord                  [m]
# lv = 10.2                       # VT tail arm                   [m]
# Av_eff = 3.5                    # VT : effective Aspect Ratio   [-]         FILE: Sizing_Vertical_Tail
# Sv = 13.96                      # Surface area VT               [m^2]
# taperv = 0.5                    # Taper VT                      [-]
# Vh_V = 1                        # Volume fraction for T-tail    [-]
# Av = 4                          # VT aspect ratio               [-]
# Sh = 11.376                     # Surface area HT               [m^2]
# taperh = 0.75                   # Taper HT                      [-]
# x_h = 23.3                      # x-location horizontal tail    [m]
# A_h = 8                         # HT aspect ratio               [-]
# lh = 10.6                       # Tail arm                      [m]