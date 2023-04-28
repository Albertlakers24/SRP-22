''' AIRCRAFT STABILITY & CONTROL OUTPUTS --> INES  ---> Last Update - 27/04/23 '''
import numpy as np

'''TO BE UPDATED:'''
AC_loc = 13                             # ac location           [m]
Cybeta_v = -0.630020450020227           # In y direction        [-]             FILE:
x_AC = 13                               # AC location           [m] aft of cg
Cnbeta = 0.09237746426846452            # Cnbeta                [-]             FILE: Vertical_Tail
Cn_r = -0.18                            # Cn_r                  [-]             FILE: Vertical_Tail
Cydelta_r = 0.25                        # Cydelta_r             [-]             FILE: Vertical_Tail
Cndelta_r = -0.09                       # Cndelta_r             [-]             FILE: Vertical_Tail
Clbeta = -0.031725211066274815          # Clbeta                [-]             FILE: Static_Derivatives
Cn_delta_a = 1                          # Cn_delta_a            [rad^-1]        FILE: ??
Cl_delta_a =1                           # Cl_delta_a            [rad^-1]        FILE: ??
CNh_delta = 0.024*180/np.pi             # Normal force gradient [rad^-1]        FILE: Sizing_Horizontal_Tail2
Chdelta = -0.020525670157042666         # Ch_delta              [rad^-1]        FILE: Sizing_Horizontal_Tail2
Chalpha = -0.017860003939365404         # Ch_alpha              [rad^-1]        FILE: Sizing_Horizontal_Tail2
CNhalpha_free = 1.3534851890535113      # Cnhalpha_free         [rad^-1]        FILE: Sizing_Horizontal_Tail2

'''INPUTS FOR THE SYMMETRIC DYNAMIC MOTIONS'''
CZadot =1                               # CZadot                [
Cmadot =1
KY2 =1.25 * 1.114
CXu= 1
CXa =1
CZ0=1
CXq=1
CZu=1
CZa=1
CX0 =1
CZq=1
Cmu=1
Cmalpha=-0.8866909436642294             # Cm_alpha              [rad^-1]        FILE: Derivative_Calculations_Dynamic
Cmq=1
CXde=1
CZde=1
Cmde =1
