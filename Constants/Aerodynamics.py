''' AERODYNAMICS OUTPUTS --> MEGHA  ---> Last Update - 08/03/22 '''

'''To Be Added:'''
R_lfus = 115629986.920      # Fuselage Reynolds number              [-]
CL_Alpha_VT =  4.03         # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT
CL_Alpha_VT_eff = 3.77      # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT
Cl_Alpha_VT_Airfoil = 1     # VT :  cl over alpha airfoil of VT     [rad^-1]
downwash = 0.00571439       # Downwash gradient                     [-]             FILE:
CL_Alpha_HT =  2.55         # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT
Cl_Alpha_HT_Airfoil = 1     # VT :  cl over alpha airfoil of VT     [rad^-1]
CLH_adj = -0.8                  # HT: CL for an adjustable tail         [-]
CD0_tailhCR = 0.0741993961878728 # HT: CD0 for cruise condition     [-]
CLwf = 0.8                  # Cruise: wing and fuselage lift coefficient    [-]

''' Zero-Drag estimations '''
CD0_CR =  0.021106443565976893      # CD0 Clean/Cruise [-]
CD0_15 = 0.04030446169943226        # CD0 Flap deflection 15 deg [-]
CD0_40 = 0.07520760765686323        # CD0 Flap deflection 40 def [-]

''' Cruise data '''
CL_DesCruise = 0.63369832                 # CL Design Cruise [-]
CD_DesCruise = 0.03323038                # CD Design Cruise [-]
CL_CD_DesCruise = 18.31272406 #19.06985102             # CL/CD Design Cruise [-]
Alpha_DesCruise = 4                       # Angle of attack Cruise [deg]
CL_Max_Clean =  1.3499704142011835        # CL max wing clean configuration [-]

'''Take-off data --> 40 deg flap deflection'''
CL_DesTakeOff = 2.1180312426              # CL max 40 deg deflection/1.13 [-]
CD_DesTakeOff = 0.17041399                # CD at 40 deg flap deflection [-]
CL_CD_TakeOff = 12.428740               # CL/CD at 40 deg flap deflection [-]
Alpha_DesTakeOff = 11.0923                # Alpha at 40 deg flap deflection [-]
CL_MaxTakeOff = 1.9                       # CL max wing at flap deflection 40 deg

''' Landing data --> 40 deg flap deflection'''
CL_DesLand = 2.1180312426              # CL max 40 deg deflection/1.13 [-]
CD_DesLand = 0.17041399                # CD at 40 deg flap deflection [-]
CL_CD_Land = 12.428740                 # CL/CD at 40 deg flap deflection [-]
Alpha_DesLand = 11.0923                # Alpha at 40 deg flap deflection [-]
CL_MaxLand = 2.2                       # CL max wing at flap deflection 40 deg

''' Slopes '''
CL_Alpha_Wing = 5.63403713             # Lift slope wing [1/rad]
Cl_Alpha_WingAirfoil =  5.99001331     # Airfoil lift slope [1/rad]
Cl_Alpha_TailAirfoil =  5.7295779513   # Tail Airfoil lift slope [1/rad]

''' TO BE UPDATED '''
DeltaCLflaps = 0.94                    # Change in CL wing due to flaps      [-]    todo: add
DeltaClmax = 2.3                       # Change in airfoil Cl by flaps       [-]    todo: add
CL0_Land = 0.29971                     # CL0 for flapped wing (landing)      [-]
Cm0_airfoil = -0.0794                  # Cm for airfoil at alpha=0           [-]