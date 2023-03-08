import numpy as np
import matplotlib.pyplot as plt
from Constants.FlightPerformance_Propulsion import eta_prop, eta_EM, eta_wire, eta_fuelcell, eta_inverter, e_lh2
from Constants.MissionInputs import h_cruise, t_loiter, V_cruise, f_con, g, R_div, R_norm
from Constants.Aerodynamics import CL_CD_DesCruise
from Constants.Masses_Locations import m_f, m_oem, m_mto, m_pldes, reserve_fuel, trip_fuel
# Propulsion charecteristics

#LD_crs = 16.787 #CL_CD_Des_CR
#eta_eng = 0.6 * 0.97 * 0.995 * 0.95

eta_eng = eta_fuelcell * eta_wire * eta_inverter * eta_EM
pl_increase = 1.04
m_plmax = m_pldes * pl_increase

def max_payload_mass(m_pldes, pl_increase, m_mto, m_oem):
    m_plmax = min(m_pldes * pl_increase, m_mto - m_oem)
    return m_plmax

def fuelmass_maxpl(m_mto, m_oem, m_plmax):
    m_fB = max(0, m_mto - m_oem - m_plmax)
    return m_fB

def R_cruise(m_pl, m_fuel):
    r_tot =  eta_eng * eta_prop * (CL_CD_DesCruise) * (e_lh2/g) * np.log((m_oem + m_pl + m_fuel)/(m_oem + m_pl))
    return r_tot

def R_tot(R_norm, R_cruise, R_div, t_loiter):
    R_lost = (1 / 0.7) * (CL_CD_DesCruise) * (h_cruise + ((V_cruise ** 2) / (2 * g)))
    R_eq = ((R_norm + R_lost) * (1 + f_con)) + (1.2 * R_div) + (t_loiter * V_cruise)
    R_aux = R_eq - R_norm
    R = R_cruise - R_aux
    return R

def plotting(ranges, plmasses, masses, colour1, colour2):
    # plotting the points
    plt.plot(ranges, plmasses, color=colour1, linewidth=3,
             marker='o', markerfacecolor=colour1, markersize=5, label = 'Payload')
    plt.plot(ranges, masses, color = colour2, linewidth=3,
             marker='o', markerfacecolor=colour2, markersize=5, label = 'Aircraft mass')

    plt.axhline(y=plmasses[1], color='grey', linestyle='--')
    plt.annotate('Design payload', xy=(1600, plmasses[1] + 60))

    plt.axvline(x=ranges[1], color='grey', linestyle='--')
    plt.annotate('Range @ Design payload', xy=(ranges[1] -60, 7000), rotation='vertical')

    plt.axvline(x=ranges[2], color='grey', linestyle='--')
    plt.annotate('Divergence Range', xy=(ranges[2] -60, 7000), rotation='vertical')

    plt.axvline(x=ranges[3], color='grey', linestyle='--')
    plt.annotate('Ferry Range', xy=(ranges[3] - 60, 7000), rotation='vertical')

    plt.axhline(y=masses[1], color='grey', linestyle='--')
    plt.annotate('Maximum take off mass', xy=(1600, masses[1] + 60))

    plt.axhline(y = m_oem, color='grey', linestyle = '--')
    plt.annotate('Operational empty mass', xy=(1500, m_oem + 50))

    plt.xlim(0,2900)
    plt.ylim(0,19500)
    n = ['A', 'B', 'C', 'D']
    for i, txt in enumerate(n):
        plt.annotate(txt, (ranges[i], plmasses[i]))
    for i, txt in enumerate(n):
        plt.annotate(txt, (ranges[i], masses[i]))
    # naming the x axis
    plt.xlabel('Range [nmi]')
    # naming the y axis
    plt.ylabel('Mass [kg]')
    # giving a title to my graph
    #plt.title('Payload range')
    plt.legend()
    # function to show the plot
    plt.show()

# -------------------- POINT A ----------------

massA = m_oem + m_pldes
plmassA = m_pldes
rangeA = 0
reserveA = plmassA+reserve_fuel
# -------------------- POINT B ----------------

massB = m_mto
plmassB = m_pldes
fuelB = trip_fuel
r_B = R_cruise(m_pldes, fuelB)
rangeB = (R_tot(r_B, r_B, 0, 0))/1852
reserveB = plmassB+reserve_fuel

# -------------------- POINT C ----------------
R_nom = R_norm*1852
massC = m_mto
plmassC = m_pldes
fuelC = reserve_fuel + trip_fuel
r_C = R_cruise(m_pldes, fuelC)
rangeC = (R_tot(r_C, r_C, R_div,t_loiter))/1852
reserveC = plmassC+reserve_fuel

# -------------------- POINT D ----------------
massD = m_oem + reserve_fuel + trip_fuel
plmassD = 0
r_D = R_cruise(0, m_f)
rangeD = (R_tot(r_D, r_D, 0, 0))/1852
reserveD = plmassD+reserve_fuel

ranges= [rangeA,rangeB, rangeC, rangeD]
masses = [massA, massB, massC, massD]
plmasses = [plmassA, plmassB, plmassC, plmassD]

plotting(ranges, plmasses, masses, 'indianred', 'forestgreen')

