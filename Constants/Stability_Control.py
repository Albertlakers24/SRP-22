''' AIRCRAFT STABILITY & CONTROL OUTPUTS --> INES  ---> Last Update - 27/04/23 '''
import numpy as np

'''TO BE UPDATED:'''
RFCO2_ATR =1.015*10**-11                    # RF for CO2            [W/m^2/s]       FILE: CO2 (graph)
RFNOx = -5.5054104*10**-16              # RF for NOx            [W/m^2/s]       FILE: NOx
AC_loc = 13                             # ac location           [m]
x_AC = 13                               # AC location           [m] aft of cg
KY2 =1.25 * 1.114                       # KY^2                  [-]             FILE: Derivative_Calculations_Dynamic
KX2= 0.02                               # KX^2                  [-]             FILE: Derivative_Calculations_Dynamic
KZ2 = 0.042                             # KZ^2                  [-]             FILE: Derivative_Calculations_Dynamic

'''UPDATED :)'''
CNh_delta = 0.024*180/np.pi             # Normal force gradient [rad^-1]        FILE: Sizing_Horizontal_Tail2
Chdelta = -0.07305017255085192          # Ch_delta              [rad^-1]        FILE: Sizing_Horizontal_Tail2
Chalpha = -0.022660754006932418         # Ch_alpha              [rad^-1]        FILE: Sizing_Horizontal_Tail2
CNhalpha_free = 3.9031332654360016      # Cnhalpha_free         [rad^-1]        FILE: Sizing_Horizontal_Tail2
Cmdelta = -1.143973535853557            # Cmdelta               [rad^-1]        FILE: Sizing_Horizontal_Tail2
Cnbeta = 0.007275695612859956           # Cnbeta                [-]             FILE: Sizing_Vertical_Tail
Cn_r = -0.1759679981035936              # Cn_r                  [-]             FILE: Sizing_Vertical_Tail
Cydelta_r = 0.46681006064945546         # Cydelta_r             [rad^-1]        FILE: Sizing_Vertical_Tail
Cndelta_r = -0.06679567439999999        # Cndelta_r             [rad^-1]        FILE: Sizing_Vertical_Tail
Cybeta_v = -0.4718537544771005          # In y direction        [-]             FILE: Sizing_Vertical_Tail
Clbeta = -0.17721899797340146           # Clbeta                [-]             FILE: Derivatives_Calculations_Static
Clr = 0.2567204919319156                # Clr                   [-]             FILE: Derivatives_Calculations_Static
Cn_delta_a = 0.12602444600489324        # Cn_delta_a            [rad^-1]        FILE: Derivatives_Calculations_Static
Cl_delta_a = -1.9887135885872829        # Cl_delta_a            [rad^-1]        FILE: Derivatives_Calculations_Static
Cybeta = -0.9511510304734292            # Cy_beta               [-]             FILE: Derivatives_Calculations_Static
Clp= -0.5012864466894098                # Clp                   [-]             FILE: Derivatives_Calculations_Static
Cnp = -0.08227687637026661              # Cnp                   [-]             FILE: Derivatives_Calculations_Static
CX0 =0                                  # CX0                   [-]             FILE: Derivative_Calculations_Dynamic
CXu= 0                                  # CXu                   [-]             FILE: Derivative_Calculations_Dynamic
CXa =0.3848768750994293                 # CXa                   [-]             FILE: Derivative_Calculations_Dynamic
CXq=0                                   # CXq                   [-]             FILE: Derivative_Calculations_Dynamic
CZ0= -0.23009078901223706               # CZ0                   [-]             FILE: Derivative_Calculations_Dynamic
CZu= -1.26739664                        # CZu                   [-]             FILE: Derivative_Calculations_Dynamic
CZa= -5.63403713                        # CZa                   [-]             FILE: Derivative_Calculations_Dynamic
CZadot = -0.8928919690218844            # CZadot                [-]             FILE: Derivative_Calculations_Dynamic
CZq= -7.203936980288712                 # CZq                   [-]             FILE: Derivative_Calculations_Dynamic
Cmu= 0                                  # Cmu                   [-]             FILE: Derivative_Calculations_Dynamic
Cmadot = -3.2328756039172455            # Cmadot                [-]             FILE: Derivative_Calculations_Dynamic
Cmq= -14.345730623699906                # Cmq                   [-]             FILE: Derivative_Calculations_Dynamic
Cmalpha= -1.5970697978937585            # Cm_alpha              [rad^-1]        FILE: Derivative_Calculations_Dynamic
muc= 96.05875842024848                  # muc                   [-]             FILE: Derivative_Calculations_Dynamic
mub= 8.388801600711354                  # mub                   [-]             FILE: Derivative_Calculations_Dynamic
CXde= 0                                 # CXdelta               [-]             FILE: Derivative_Calculations_Dynamic
CZde= -0.3159554861002181               # CZdelta               [-]             FILE: Derivative_Calculations_Dynamic