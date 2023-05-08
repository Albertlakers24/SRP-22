import numpy as np
import matplotlib.pyplot as plt
from Constants.Stability_Control import RFNOx, RFCO2_ATR

"""TODO after the holidays:
* add massflow as a function to other files
* check for fleet size with marketing group
* check if forcing factor should be made a function with time? -> Schmidt-Appleman criterion
* check RF parameter for contrails if changed in case of a hydrogen Fuel Cell
"""


"""CONSTANTS FROM OTHER PARAMTERS"""
"""horizon is the time horizon- 30yrs (Dallara) for short-lived species and 100 yrs for long-lived species SRP-22 and ATR"""
horizon = 30                            # Time horizon                              [years]
massflow = 0.041                        # Average massflow engine                   [kg/s]
t_year = 3098*60*60                     # Time flown by 1 aircraft for 1 year       [sec]
number_aircraft = 25*1000/30            # Number of aircraft 1,000 in 30 years with a lifetime of 25 years     [number of aircraft]
Distance_year = 703209                  # Distance flown by 1 aircraft for 1 year   [nmi]

massflow_ATR = 143.89*10.84116*0.549*0.25 # Average massflow engine                   [kg/s]

"""CONSTANT PARAMTERS """
RF_CO2 = 3.7                            # Radiative Forcing CO2, doubling effect    [W/m^2]
S = 3                                   # Climate sensitivity parameter             [K]
alpha_t = 0.595                         # [-]
tau_t1 = 8.4                            # Time constant                             [yrs]
tau_t2 = 409.5                          # Time constant                             [yrs]
Eff_contrails = 0.59                    # Efficacy for contrails                    [-]
forcing_factor_contrails = 0.83         # Forcing factor for contrails              [-]
RF_ref_over_Lref = 2.21*10**-12         # RF_ref/L_ref                              [(W/m^2)/nmi]
EI_water = 2.6 * 1.26                   # Emission Index                            [kg/kg]
ratio_water = 7.43 * 10 ** -15          # RFref/Eref for water                      []
Eff_water = 1.14                        # Efficacy for water                        [-]

"""TIMES"""
time = np.arange(2035, 2035+horizon, 1)
time_horizon = np.arange(0,horizon, 1)


def convolve(signal, response):
    """ Source: Proesmans code
    Determines the convolution integral for the given signal and response (or
    filter) function.
    :param signal: Signal, such as emissions over years
    :param response: Response function to be applied to signal, such as decay
    or other filter
    :return: Convolution integral of provided signal and response function
    """
    output = np.convolve(signal, response, 'full')
    return output[:signal.size]

def RF_contrails(RF_CO2):
    """ Source: Proesmans
    :param lenght: dependent on mission  50,000 (nmi/year)          -DONE
    :param forcing_factor_contrails: dependent on altitude (-)      -TO BE REWRITTEN DEPENDENT ON ALTITUDE
    :return: RF_contrail_norm: normalized RF (W/m^2/year)           -DONE
    """
    length = Distance_year*number_aircraft                          # input: [nmi/year] for a fleet
    t_test = np.ones(horizon)
    RF_contrail = RF_ref_over_Lref*length*forcing_factor_contrails
    RF_contrail_norm = t_test*Eff_contrails * RF_contrail / RF_CO2
    return RF_contrail_norm

def RF_water(massflow, t_year, RF_CO2):
    """ Source: Proesmans
    Calculation of the normalized RF for water emissions
    :param ratio: fixed RF parameter for water (W/m^2/kg)           -DONE
    :param EI: Emission index for water (kg/kg)                     -DONE
    :param massflow: engine (kg/s)                                  -DONE
    :param RF_CO2: doubling coefficient (W/m^2)                     -DONE
    :return: normalized RF (W/m^2/year)                             -DONE
    """
    mass = massflow*t_year*number_aircraft                  # output: [kg/year]
    E = EI_water*mass                                       # input: [kg/year] for a fleet
    t_test = np.ones(horizon)
    RF_water = ratio_water*E
    RF_water_norm = t_test*Eff_water*RF_water/RF_CO2
    return RF_water_norm

t_test = np.ones(horizon)
print("t_test",t_test)

def G_factor(time):
    """ Source: https://arc.aiaa.org/doi/pdf/10.2514/1.J050763
    :param time: time array (years)                                   -DONE
    :return: G-factor (/year)                                         -DONE
    """
    G_factor = S*((alpha_t/tau_t1)*np.exp(-time/tau_t1) + ((1-alpha_t)/tau_t2)*np.exp(-time/tau_t2))
    return G_factor

def Temperature_response(RF, time):
    """
    Estimates of temperature response
    :param RF is the total RF for contrails and water
    :param time is the time in years
    :return: Temperature response (K/year)
    """
    result = convolve(signal=RF,response=G_factor(time=time_horizon))
    return result

"""RF_total is dependent on time (per years)"""
RF_total = RF_contrails(RF_CO2=RF_CO2) + RF_water(massflow=massflow, t_year=t_year, RF_CO2=RF_CO2)
print("RF_total=", RF_total, "per year")
print("RF_con=", RF_contrails(RF_CO2=RF_CO2))
print("RF_water=", RF_water(massflow=massflow, t_year=t_year, RF_CO2=RF_CO2))

"""Result on Delta T"""
result_Temperatureresponse = Temperature_response(RF=RF_total,time=time_horizon)
print("Temperature Change=",result_Temperatureresponse)

"""For ATR:"""
RFCO2_ATR = (1/0.6931)*np.log((380+0.4)/380)
RF_total_ATR = RF_contrails(RF_CO2)+RF_water(massflow,t_year,RF_CO2)+RFNOx+RFCO2_ATR
result_Temperatureresponse_ATR = Temperature_response(RF=RF_total_ATR,time = time_horizon)
print("RF_total_ATR =", RF_total_ATR)
print("RF_CO2=",RFCO2_ATR)

fig, ax = plt.subplots()

# plt.plot(time, result_Temperatureresponse)
ax.semilogy(time, result_Temperatureresponse,label="SRP-22")
ax.semilogy(time, result_Temperatureresponse_ATR,label="ATR")
plt.grid(True)
plt.xlabel("Time [year]", fontsize=15, weight="bold")
plt.ylabel("Average Temperature Change [K]", fontsize=15, weight="bold")
plt.legend(loc="lower right")
# plt.title("Average Temperature Response for a time horizon of 30 years for an entire fleet", fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

# plt.plot(time, result_Temperatureresponse_ATR)
# plt.grid(True)
# plt.xlabel("Time [year]", fontsize=15, weight="bold")
# plt.ylabel("Average Temperature Change [K]", fontsize=15, weight="bold")
# # plt.title("Average Temperature Response for a time horizon of 30 years for an entire fleet", fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.show()

# plt.plot(time, G_factor(time=time_horizon))
# plt.grid(True)
# plt.xlabel("Year")
# plt.ylabel("Scaled temperature impulse response")
# plt.title("Scaled Temperature Impulse Response function for designated time")
# plt.show()

#RF_test needs to be made an array
print("length RF=", len(RF_total))
print("length time", len(time))
print("length G=", len(G_factor(time=time)))
print("type RF=", type(RF_total))
print("type time=", type(time))
print("type G_factor=", type(G_factor(time=time)))

