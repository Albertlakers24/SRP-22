''' FLIGHT PERFORMANCE AND PROPULSION OUTPUTS --> BEREND  ---> Last Update - 27/2/22 '''

# Suggested headers - pls change if needed -- segregate as you'd like
''' FLIGHT PERFORMANCE '''
n_ult_pos = 3.24796708 * 1.5
n_ult_neg = -1.24796708 * 1.5

''' PROPULSION SYSTEM '''

''' Power '''
Power_tot = 2473.981955 * 10**3        # Total required power [W]
Voltage_cell = 0.7              # Cell volatge [V]

''' Hydrogen tank '''


''' Efficiencies '''
eta_prop = 0.85                     #Propeller efficiency
eta_EM = 0.95                       #Electric motor efficiency
eta_wire = 0.97                     #Wire efficiency
eta_inverter = 0.995                #Inverter efficiency
eta_fuelcell = 0.60                 #Fuel cell efficiency
fc_power_density = 3                #kW/kg - fuel cell power density
inverter_power_density = 30         #kW/kg - inverter power density
em_power_density = 8                #kW/kg - electric motor power density
H2_power_density = 120              #MJ/kg
e_lh2 = 120*10**6                   #J/kg Specific Energy Liquid Hydrogen
