import numpy as np
import pandas as pd
import os
from Constants.AircraftGeometry import S_w,taperw,y_mac
from Aerodynamics.HLD_choices import HLD_TE_deltaClmax,TE_hinge_line_angle_deg
from Constants.Aerodynamics import CL_Alpha_Wing

working_directory = os.getcwd()
working_directory = working_directory.split("SRX-22-ALAR-FOX")[0]
working_directory = os.path.join(working_directory,"SRX-22-ALAR-FOX")
folder = os.path.join(working_directory,"Aerodynamic_characteristics")
file = os.path.join(folder,"HLD_design_choices.txt")
HLD_choices = pd.read_csv(file,sep = "\t")
HLD_choices = HLD_choices[HLD_choices['ld[m]']>0].sort_values("ld[%]")
single_slotted_HLD_choices = HLD_choices[HLD_choices['flap_type']== "single slotted+None"]
double_slotted_HLD_choices = HLD_choices[HLD_choices['flap_type']== "double slotted+None"]
fowler_HLD_choices = HLD_choices[(HLD_choices['df'] < 43) & (HLD_choices["flap_type"] == "fowler+None") & (HLD_choices["Cf_C"]>=0.35) & (HLD_choices["df"]==40)]
HLD_final_design = HLD_choices.iloc[2]
lm = 1.5
ld = HLD_final_design[0]
ld_perc = HLD_final_design[1]
SwfS = HLD_final_design[2]
Cf_C = HLD_final_design[3]
DCL = HLD_final_design[4]
DPA = HLD_final_design[5]
DCl = HLD_final_design[6]
flap_deflection = HLD_final_design[7]
c_prime_over_c_land = HLD_final_design[8]
S_total_land = S_w * (1 + SwfS*(c_prime_over_c_land-1))
cf_over_cprime = Cf_C * (1/c_prime_over_c_land)
c_prime_over_c_to = (0.4 + 15 * (0.2/(45-(1/1.2*10))))*Cf_C + 1
S_total_to = S_w * (1 + SwfS*(c_prime_over_c_to-1))
DCL_to = 0.9 * HLD_TE_deltaClmax(Cf_C,15,"double slotted")[0] * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
DPA_to = -10 * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg))
DCL_alpha_land = S_total_land/S_w * CL_Alpha_Wing
DCL_alpha_to = S_total_to/S_w * CL_Alpha_Wing
# print(DCL_to)
# print(DPA_to)
# print(DCL_alpha_to)
