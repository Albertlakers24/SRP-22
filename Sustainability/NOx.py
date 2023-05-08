import numpy as np


def convolve(signal, response):
    """
    Determines the convolution integral for the given signal and response (or
    filter) function.
    :param signal: Signal, such as emissions over years
    :param response: Response function to be applied to signal, such as decay
    or other filter
    :return: Convolution integral of provided signal and response function
    """

    output = np.convolve(signal, response, 'full')
    return output[:signal.size]


def concentration_response(t_test,
                           alpha_i=np.array(
                               [0.067, 0.1135, 0.152, 0.097, 0.041]),
                           tau_i=np.array([np.inf, 318.8, 79.8, 18.8, 1.7])):
    """
    Delta CO2 concentration response due to CO2 emissions.
    :param t_years: Array with time, typically years
    :type t_years: np.ndarray
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: Delta CO2 concentration response function
    """
    tn = t_test - t_test[0]  # + 1
    alpha_in = alpha_i * 1e-3 * 1e-9  # convert from ppbv/Tg to ppmv/kg
    try:
        one = np.ones((5, tn.size)) * -1.
    except AttributeError:
        one = np.ones((5, 1)) * -1.
    one = one * tn
    one = one.transpose()
    one = one / tau_i
    one = np.exp(one)

    one = one * alpha_in

    one = np.sum(one, axis=1)
    return one


#############################################################################
#For NOx emissions
#Assumption for s=1

t_test = np.arange(0,8974800,1)          #s - 1 year  - 4x500 missions
Type = 1
#1-NOx -> CH4
#2-NOx -> O3L
Configuration = 4
#0 - kerosene
#1 - Hydrogen Combustion

tau_i = 12*365*24*60    #sec

if Type == 1:
    A_i = -5.16 * 10 ** -13  # W/m2/kg
if Type == 2:
    A_i = -1.21 * 10 ** -13  # W/m2/kg
if Configuration ==0:
    massflow=0.1805555            #kg/s
    EI_NOx = 0.014          #kg/kg
if Configuration==1:    #Hydrogen combustion
    massflow = 0.05843          # kg/s
    EI_NOx = 0.0028         # kg/kg
if Configuration==3:    #HES
    massflow = 0.092525739          # kg/s
    EI_NOx = 0.014          # kg/kg
if Configuration==4:    #EHSP
    massflow = 0.07621          # kg/s
    EI_NOx = 0.014          # kg/kg


def emissions_NOx(EI_NOx, massflow):
    """
    :param EI_NOx: emission index (kg/kg)
    :param massflow: mass flow of respective configuration (kg/s)
    :return: emissions NOx (kg/s)
    """
    output_NOx = EI_NOx*massflow
    return output_NOx

emissions_NOx(EI_NOx=EI_NOx, massflow=massflow)
#print("emissions NOx=", emissions_NOx(EI_NOx=EI_NOx, massflow=massflow, t_year=t_test))
emissionsNOx = np.sum(emissions_NOx(EI_NOx=EI_NOx, massflow=massflow))

def Gfactor (t_year, A_i, tau_i):
    """
    :param t_year: 1 year 4x500 nmi missions (s)
    :param A_i: W/m2/kg
    :param tau_i: in sec
    :return:W/m2/kg
    """
    G = A_i*np.exp(-t_year/tau_i)
    return G

#print("Gfactor=",Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i))

def concentration_NOx(emissionsNOx, t_test, A_i,tau_i):
    """
    Estimates CO2 concentration in atmosphere in ppmv.
    :param t_year: Array with time, typically years
    :type t_year: np.ndarray
    :param emissions_co2: CO2 emissions per year in t_year, in kg
    :type emissions_co2: np.ndarray
    :param co2_0: Pre-industrial or background CO2 concentration in
        ppmv
    :type co2_0: float
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: concentration NOx in
    """
    result = convolve(signal=emissionsNOx,response=Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i))
    return result

#Print statements
print("Sum of NOx emissions=",emissionsNOx,"kg")
print("Sum of Gfactor=", np.sum(Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i)))
print("RFNOx =", np.sum(concentration_NOx(emissionsNOx=emissionsNOx, t_test=t_test, A_i=A_i, tau_i=tau_i)), "W/m^2/s")
RFNOx = np.sum(concentration_NOx(emissionsNOx=emissionsNOx, t_test=t_test, A_i=A_i, tau_i=tau_i))