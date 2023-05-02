''' AERODYNAMICS OUTPUTS --> MEGHA  ---> Last Update - 02/05/22 '''

'''To Be Added:'''
R_HT =7.5*10**6                     # Reynolds number of the HT
Cmac = -0.6926729205077256          # Cmac                          [-]             FILE: ScissorPlot_HT
R_lfus = 111301293.64412345  # Fuselage Reynolds number              [-]
CL_Alpha_VT =  2.102#4.03         # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT TODO: Check
CL_Alpha_VT_eff = 3.77      # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT TODO: Check
Cl_Alpha_VT_Airfoil = 5.841193106203845      # VT :  cl over alpha airfoil of VT     [rad^-1]
downwash = 0.24789#0.00571439       # Downwash gradient                     [-]             FILE: TODO: Update/Check
CL_Alpha_HT =  4.3297#2.55         # VT : CL over alpha                    [rad^-1]        FILE: ScissorPlot_HT  TODO: Check
Cl_Alpha_HT_Airfoil = 5.841193106203845      # VT :  cl over alpha airfoil of VT     [rad^-1]
CLH_adj = -0.8                  # HT: CL for an adjustable tail         [-]
CD0_tailhCR = 0.09292522985159975 # HT: CD0 for cruise condition     [-]
CLwf = 0.63369832                  # Cruise: wing and fuselage lift coefficient - steady state lift coeff at 1g flight   [-]

''' Zero-Drag estimations '''
CD0_CR = 0.0227181      # CD0 Clean/Cruise [-]
CD0_15 = 0.042574742757        # CD0 Flap deflection 15 deg [-]
CD0_40 = 0.017556784573        # CD0 Flap deflection 40 def [-]

''' Cruise data '''
CL_DesCruise = 0.63369832                 # CL Design Cruise [-]
CL_CD_DesCruise = 18.28261      # CL/CD Design Cruise [-]
Alpha_DesCruise = 3.962389                # Angle of attack Cruise [deg]
CL_Max_Clean =  1.3499704142011835        # CL max wing clean configuration [-]
CD_DesCruise = 0.034661
'''Take-off data --> 15 deg flap deflection'''
CL_DesTakeOff = 1.9              # CL max 40 deg deflection/1.13 [-]
CD_DesTakeOff = 0.1331683                # CD at 40 deg flap deflection [-]
CL_CD_TakeOff = 13.068810               # CL/CD at 40 deg flap deflection [-]
CL_MaxTakeOff = 1.9                       # CL max wing at flap deflection 40 deg

''' Landing data --> 40 deg flap deflection'''
CL_DesLand = 2.1180312426              # CL max 40 deg deflection/1.13 [-]
CD_DesLand = 0.19790                # CD at 40 deg flap deflection [-]
CL_CD_Land = 10.6112                 # CL/CD at 40 deg flap deflection [-]
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