import numpy as np
import matplotlib.pyplot as plt
from Constants.Masses_Locations import m_mto,W_S_design
from Constants.MissionInputs import ft_m, nmi_m, kts_m_s,ISA_calculator, rho_0,V_cruise,rho_cruise,g, rho_5000
from Constants.Aerodynamics import CL_MaxTakeOff, CL_CD_TakeOff, CD0_15, CL_CD_DesCruise, CD_DesCruise, CL_DesCruise, Oswald_TO, CL_MaxLand
from Constants.FlightPerformance_Propulsion import H2_power_density
from Constants.AircraftGeometry import S_w, Aw
import sympy as sp
import math

e_take_off = Oswald_TO
C_LFL = 0.45
m_mto = m_mto
m_fuel = 611
h2 = 50 * ft_m
P_available_prop = 2.5 * 10**6
s_0 = 0
s_1 = 5000 * ft_m
s_2 = 15000 * ft_m
s_3 = 28000 * ft_m
s_4 = 10000 * ft_m
s_5 = 0
ROC1 = 6
ROC2 = 5
ROC3 = 5
ROD1 = -4.7386#-1500 * ft_m / 60
ROD2 = -3.77#-1110 * ft_m / 60

V_IAS1 = 140 * kts_m_s
V_IAS2 = 176 * kts_m_s
V_IAS3 = 176 * kts_m_s
V_IAS4 = 176 * kts_m_s
V_IAS5 = 140 * kts_m_s
V_IAS_list = [V_IAS1, V_IAS2, V_IAS3]
s_list = [s_0, s_1, s_2, s_3, s_4, s_5]
V_TAS_list = []
def V_TAS_calc(V_IAS, h, h_0):
    # V_TAS = V_IAS * (1 + (h + h_0)/ 2 / 1000 * 0.02)
    avg_h = (h + h_0) / 2
    T, p, rho, a = ISA_calculator(avg_h, 0)
    V_TAS = V_IAS * np.sqrt(rho_0 / rho)
    V_TAS_list.append(V_TAS)
    return V_TAS

V_TAS1_avg = V_TAS_calc(V_IAS1, s_1, s_0)
V_TAS2_avg = V_TAS_calc(V_IAS2, s_2, s_1)
V_TAS3_avg = V_TAS_calc(V_IAS3, s_3, s_2)
V_TAS4_avg = V_TAS_calc(V_IAS4, s_3, s_4)
V_TAS5_avg = V_TAS_calc(V_IAS5, s_4, s_5)

ROC_list = [ROC1, ROC2, ROC3, ROD1, ROD2]

V_x = []
for i in np.arange(0, len(V_TAS_list), 1):
    V_x_part = np.sqrt(V_TAS_list[i]**2 - ROC_list[i]**2)
    V_x.append(V_x_part)

def t_calc():
    t_list = []
    for i in np.arange(0, 3):
        t = (s_list[i+1] - s_list[i]) / ROC_list[i]
        t_list.append(t)
    for i in np.arange(3,5):
        t = (s_list[i] - s_list[i+1]) / ROC_list[i] * -1
        t_list.append(t)
    return t_list
t_list = t_calc()
descent_distance = V_x[3]*t_list[3]+V_x[4]*t_list[4]
total_distance = sum(np.array(V_x) * np.array(t_list))
full_time_climb = t_list[0] + t_list[1] + t_list[2]
full_time_descent = t_list[3] + t_list[4]
print(f"Full time of descent: {full_time_descent / 60} minutes")
print(f"Descent distance: {descent_distance / nmi_m} nmi")
cruise_distance = (1000 * 1852 - total_distance)                             #in m
cruise_time = cruise_distance / V_cruise
full_time = full_time_climb + cruise_time + full_time_descent

Weight = m_mto
C_L = 0.72
Lift = 1/2 * rho_cruise * S_w * C_L * V_cruise**2
CL_need = (Weight * 9.80665) / (1/2 * rho_cruise * S_w * V_cruise**2)
S_opt = (Weight * 9.80665) / (1/2 * rho_cruise * C_L * V_cruise**2)

# THRUST AND POWER CALCULATIONS FOR PROPELLERS
# S_big_prop = np.pi / 4 * big_prop_diameter ** 2
CL_here = 0.6296
CD_here = CL_here / 16.787

C_L_cruise = 2*0.95*W_S_design / (rho_cruise * V_cruise**2)

def mass_flow(V, h_start, h_end, dt, CL, ROC):
    h = (h_start + h_end) / 2
    T, p, rho, a = ISA_calculator(h, dt)
    CD = CL / 15
    D = 0.5 * rho * V ** 2 * S_w * CD
    Power_req = D * V
    Power_available = ROC * m_mto * g + Power_req
    mass_flow = Power_available / (H2_power_density * 10**6)
    return mass_flow, Power_available

climb1, Pa1 = mass_flow(V_TAS1_avg, s_0, s_1, 0, CL_here, ROC1)
climb2, Pa2 = mass_flow(V_TAS2_avg, s_1, s_2, 0, CL_here, ROC2)
climb3, Pa3 = mass_flow(V_TAS3_avg, s_2, s_3, 0, CL_here, ROC3)
cruise, Pa4 = mass_flow(V_cruise, s_3, s_3, 0, CL_here, 0)
descent1, Pa5 = mass_flow(V_TAS4_avg, s_3, s_4, 0, CL_here, ROD1)
descent2, Pa6 = mass_flow(V_TAS5_avg, s_4, s_5, 0, CL_here, ROD2)
Power_available = [Pa1, Pa2, Pa3, Pa4, Pa5, Pa6]
fuel_flows = [climb1, climb2, climb3, cruise, descent1, descent2]
total_fuel = 327
total_fuel_flows = sum(fuel_flows)
ff_1 = climb1 / climb1
ff_2 = climb2 / climb1
ff_3 = climb3 / climb1
ff_4 = cruise / climb1
ff_5 = descent1 / climb1
ff_6 = descent2 / climb1
ffs = [ff_1, ff_2, ff_3, ff_4, ff_5, ff_6]

t_list_all = [t_list[0], t_list[1], t_list[2], cruise_time, t_list[3], t_list[4]]
first_ff = total_fuel / sum(np.array(t_list_all) * np.array(ffs))
second_ff = ff_2 * first_ff
third_ff = ff_3 * first_ff
fourth_ff = ff_4 * first_ff
fifth_ff = ff_5 * first_ff
sixth_ff = ff_6 * first_ff
all_ffs = [first_ff, second_ff, third_ff, fourth_ff, fifth_ff, sixth_ff]
total_fuel_check = sum(np.array(t_list_all) * np.array(all_ffs))

#PEAK POWER
Peak_power = 2.5 * 10**6

#Check Mission Profile
def CL_calc(V_IAS, h, mass):
    T, p, rho, a = ISA_calculator(h, 0)
    V = V_IAS * np.sqrt(rho_0 / rho)
    CL = 2 * mass * g / (rho * V**2 * S_w)
    return CL, V, rho

def Power_available(V_IAS, h, mass, ROC):
    CL_check, V, rho = CL_calc(V_IAS, h, mass)
    C_D = CL_check / 16
    D = 1/2 * rho * V**2 * S_w * C_D
    Pa = D * V + ROC * mass * g
    if Pa < (2.5 * 10**6):
        Pa_check = "PASS"
    else:
        Pa_check = "FAIL"
    return Pa, Pa_check, CL_check

powers = []
CLs = []
i_list = []
for i in np.arange(0, 5000 * ft_m):
    Pa, P, CL_check = Power_available(V_IAS1, i, m_mto, ROC1)
    powers.append(P)
    CLs.append(CL_check)
for i in np.arange(5000 * ft_m, 15000 * ft_m):
    Pa, P, CL_check = Power_available(V_IAS2, i, m_mto, ROC2)
    powers.append(P)
    CLs.append(CL_check)
for i in np.arange(15000 * ft_m, 28000 * ft_m):
    Pa, P, CL_check = Power_available(V_IAS3, i, m_mto, ROC3)
    powers.append(P)
    CLs.append(CL_check)

Cruise_P, Cruise_Pass, CL_cruise_check = Power_available(V_cruise, 28000 * ft_m, m_mto, 0)

for i in np.arange(10000 * ft_m, 28000 * ft_m):
    Pa, P, CL_check = Power_available(V_IAS4, i, m_mto, ROD1)
    powers.append(P)
    CLs.append(CL_check)
for i in np.arange(0, 10000 * ft_m):
    Pa, P, CL_check = Power_available(V_IAS5, i, m_mto, ROD2)
    powers.append(P)
    CLs.append(CL_check)

i_list = []
for i in np.arange(0, (28000 * ft_m * 2 + 1)):
    i_list.append(i)

#PROPELLER CALCULATIONS
MTOW = m_mto * g
V_LOF = np.sqrt(2 * MTOW / (rho_5000 * S_w * CL_MaxTakeOff)) * 1.05
def Power_req_calc(V, h, dt, CL_CD):
    rho = ISA_calculator(h, dt)[2]
    CL = 2 * MTOW / (S_w * V**2 * rho)
    CD = CL / CL_CD #np.mean(find_CL_CD(CL))
    D = 1/2 * rho * V**2 * CD * S_w
    Pr = D * V
    return Pr, rho

def prop_eff(V, h, dt, CL_CD):
    eff = Power_req_calc(V, h, dt, CL_CD)[0] / P_available_prop
    return eff

def ROC_calc(V, h, dt, extra_acc, CL_CD):
    Pr, rho1 = Power_req_calc(V, h, dt, CL_CD)
    H_diff = 5
    rho2 = ISA_calculator(h+H_diff, dt)[2]
    #acceleration = V / g * np.sqrt((V * (np.sqrt(rho_0 / rho2) - np.sqrt(rho_0 / rho1)) / H_diff)**2 + extra_acc**2)
    acceleration_1 = V / g * (V * (np.sqrt(rho_0 / rho2) - np.sqrt(rho_0 / rho1)) / H_diff)
    acceleration_2 = V / g * extra_acc
    ROC = ((P_available_prop - Pr) / MTOW - acceleration_2) / (1 + acceleration_1) #/ (1 + acceleration)
    return ROC

ROC_TO = ROC_calc(57.5, 7600*ft_m, 41, 0, CL_CD_TakeOff)

def VTAS_calc(V, h, dt):
    rho = ISA_calculator(h, dt)[2]
    V_TAS = V * np.sqrt(rho_0 / rho)
    return V_TAS

V_EAS1 = 140 * kts_m_s#138.71 * kts_m_s #140
V_EAS2 = 140 * kts_m_s #176
V_EAS3 = (V_EAS2 + V_EAS1) / 2

TV2 = sp.Symbol("TV2")
CL2 = CL_MaxTakeOff / 1.13**2
eq = MTOW**2 / (rho_5000 * g * S_w * CL2 * 0.85 * TV2) + 2 * h2 * (TV2 / MTOW - (CD0_15 / CL2) - CL2 / (np.pi * Aw * e_take_off)) ** (-1) - 4500 * ft_m
TV2 = sp.solve(eq, TV2)[1]
l_run = (MTOW**2 / (rho_5000 * g * S_w * CL2 * 0.85 * TV2)) * 0.8
l_air = 2 * h2 * (TV2 / MTOW - (CD0_15 / CL2) - CL2 / (np.pi * Aw * e_take_off))**(-1)
print((l_run + l_air)/ft_m, "TAKEOFF DISTANCE")
avg_a = V_LOF**2 / (2*l_run)
ROC_TO = (P_available_prop - 1 / np.mean(CL_CD_TakeOff) * MTOW) / MTOW - avg_a * V_LOF / g
print(V_LOF/avg_a+5, "takeoff time")

#HOT AND HIGH TAKEOFF
TV2_HOT = sp.Symbol("TV2_HOT")
T_hot = ISA_calculator(5000*ft_m, 41)[0]
rho_hot = ISA_calculator(5000*ft_m, 41)[2]
eq2 = MTOW**2 / (rho_hot * g * S_w * CL2 * 0.85 * TV2_HOT) + 2 * h2 * (TV2_HOT / MTOW - (CD0_15 / CL2) - CL2 / (np.pi * Aw * e_take_off)) ** (-1) - 4500 * ft_m
TV2_HOT = sp.solve(eq2, TV2_HOT)[1]
l_run_hot = (MTOW**2 / (rho_hot * g * S_w * CL2 * 0.85 * TV2_HOT))
l_air_hot = 2 * h2 * (TV2_HOT / MTOW - (CD0_15 / CL2) - CL2 / (np.pi * Aw * e_take_off))**(-1)
print((l_run_hot + l_air_hot)/ft_m, "TAKEOFF DISTANCE HOT CLIMATE")
V_LOF_HOT = np.sqrt(2 * MTOW / (rho_hot * S_w * CL_MaxTakeOff)) * 1.05
avg_a = V_LOF_HOT**2 / (2*l_run_hot)
ROC_TO = (P_available_prop - 1 / np.mean(CL_CD_TakeOff) * MTOW) / MTOW - avg_a * V_LOF / g
print(V_LOF/avg_a+5, "takeoff time HOT")

ROCS = []
heights = []
V_input = V_LOF
ROC = ROC_TO
i_step = 10
extra_acc = 0.5
V_check = []
j = 0
s_traveled = []
s_trav = 0
j_check = []
time_climb = []
CL_endurance = 0.97786
CD_endurance = 0.048984
while V_input < V_EAS1:
    V = VTAS_calc(V_input, j * ft_m, 0)
    if j == 0:
        ROC = ROC_TO
        CL_CD = CL_CD_TakeOff
    else:
        ROC = ROC_calc(V, j * ft_m, 0, extra_acc, CL_CD_DesCruise)
    time = i_step * ft_m / ROC
    Vx = math.sqrt(V**2 - ROC**2)
    Sx = Vx * time
    s_trav = s_trav + Sx
    s_traveled.append(s_trav)
    ROCS.append(ROC)
    heights.append(j)
    V_input = V_input + time * extra_acc
    V_check.append(V_input)
    time_climb.append(time)
    i = i + i_step
    j = j + i_step
for i in np.arange(j, 4010, i_step):
    if j < i < 4000:
        V = VTAS_calc(V_EAS1, i * ft_m, 0)
        ROC = ROC_calc(V, i * ft_m, 0, 0, CL_CD_DesCruise)
        time = i_step * ft_m / ROC
        Vx = math.sqrt(V ** 2 - ROC ** 2)
        Sx = Vx * time
        s_trav = s_trav + Sx
        s_traveled.append(s_trav)
        ROCS.append(ROC)
        heights.append(i)
        V_check.append(V)
        time_climb.append(time)
        j_new = 4000
if 4000 <= i < 5000:
    while V < V_EAS2:
        V_2 = VTAS_calc(V, i * ft_m, 0)
        new_acc = 0.2
        ROC = ROC_calc(V_2, i * ft_m, 0, new_acc, CL_CD_DesCruise)
        time_step = i_step * ft_m / ROC
        Vx = math.sqrt(V ** 2 - ROC ** 2)
        Sx = Vx * time
        s_trav = s_trav + Sx
        s_traveled.append(s_trav)
        ROCS.append(ROC)
        heights.append(i)
        V = V + time_step * new_acc
        i = i + i_step
        V_check.append(V)
        time_climb.append(time)
        j_new = j_new + i_step
        j_check.append(j_new)
        if i == 4990:
            break
    for i in np.arange(j_new, 5010, i_step):
        V = VTAS_calc(V_EAS2, i * ft_m, 0)
        ROC = ROC_calc(V, i * ft_m, 0, 0, CL_CD_DesCruise)
        time = i_step * ft_m / ROC
        Vx = math.sqrt(V ** 2 - ROC ** 2)
        Sx = Vx * time
        s_trav = s_trav + Sx
        s_traveled.append(s_trav)
        ROCS.append(ROC)
        heights.append(i)
        V_check.append(V)
        time_climb.append(time)
for i in np.arange(5010, 40010, i_step):#28010
    if 5010 <= i < 15000:
        V = VTAS_calc(V_EAS2, i * ft_m, 0)
        ROC = ROC_calc(V, i * ft_m, 0, 0, CL_CD_DesCruise)
        time = i_step * ft_m / ROC
        Vx = math.sqrt(V ** 2 - ROC ** 2)
        Sx = Vx * time
        s_trav = s_trav + Sx
        s_traveled.append(s_trav)
        ROCS.append(ROC)
        heights.append(i)
        V_check.append(V)
        time_climb.append(time)
    elif 15000 <= i < 28001:
        V = VTAS_calc(V_EAS2, i * ft_m, 0)
        ROC = ROC_calc(V, i * ft_m, 0, 0, CL_CD_DesCruise)
        time = i_step * ft_m / ROC
        Vx = math.sqrt(V ** 2 - ROC ** 2)
        Sx = Vx * time
        s_trav = s_trav + Sx
        s_traveled.append(s_trav)
        ROCS.append(ROC)
        heights.append(i)
        V_check.append(V)
        time_climb.append(time)

print(sum(time_climb)/60, "CLIMB TIME")

dup = {x for x in heights if heights.count(x) > 1}

#MAX SPEED CALC
V_max = ((2 * P_available_prop) / (rho_cruise * CD_DesCruise * S_w))**(1/3)
print(f"Maximum attainable velocity during cruise is {V_max/kts_m_s} kts")

#LANDING DISTANCE
s_land = m_mto * g * C_LFL * 2 / (rho_5000 * CL_MaxLand * S_w)
print(s_land/ft_m, "MAXIMUM LANDING DISTANCE")

#ENERGY CALCULATIONS
V_initial_descent = VTAS_calc(176*kts_m_s, 19000*ft_m, 0)
V_second_descent = VTAS_calc(140*kts_m_s, 5000*ft_m, 0)
glide_angle = -4.5 * (np.pi / 180)
ROD1_avg = -V_initial_descent * np.sin(glide_angle)
ROD2_avg = -V_second_descent * np.sin(glide_angle)
time_descent_1 = 18000 * ft_m / ROD1_avg
time_descent_2 = 10000 * ft_m / ROD2_avg


dist_descent = time_descent_1 * V_initial_descent * np.cos(glide_angle) + time_descent_2 * V_second_descent * np.cos(glide_angle)
total_descent_avg = time_descent_1 + time_descent_2
climb_time = sum(time_climb)
cruise_dist_avg = 1000 * nmi_m - s_traveled[2799] - dist_descent
print(cruise_dist_avg/nmi_m, "DISTANCE CRUISE")
time_normal_cruise = cruise_dist_avg / V_cruise
time_fast = cruise_dist_avg / V_max
power_cruise = 1/2 * rho_cruise * V_cruise**3 * CD_DesCruise * S_w
E_normal_cruise = power_cruise * time_normal_cruise
power_cruise_fast = 1/2 * rho_cruise * V_max**3 * CD_DesCruise * S_w
E_fast = power_cruise_fast * time_fast
print(power_cruise/10**6, "POWER NORMAL")
print(power_cruise_fast/10**6, "POWER FAST")
print(E_normal_cruise/10**6, "ENERGY NORMAL")
print(E_fast/10**6,"NORMAL FAST")

plt.plot(np.array(s_traveled) / nmi_m, heights)
plt.axhline(y=s_1 / ft_m, color='grey', linestyle='--')
plt.annotate('Start Second Climb Phase', xy=(30, s_1 / ft_m + 500))
plt.axhline(y=s_2 / ft_m, color='grey', linestyle='--')
plt.annotate('Start Third Climb Phase', xy=(2, s_2 / ft_m + 500))
plt.axhline(y=s_3 / ft_m, color='grey', linestyle='--')
plt.annotate('Start Cruise', xy=(2, s_3 / ft_m + 500))
plt.axhline(y=4000, color="grey", linestyle='--')
plt.annotate('Start acceleration no. 2', xy=(30, 4000 + 200))
plt.axhline(y=j, color="grey", linestyle="--")
plt.annotate('End acceleration no. 1', xy=(30, j + 500))
plt.axvline(x=s_trav / nmi_m, color='grey', linestyle='--')
plt.plot([s_trav / nmi_m, s_trav / nmi_m], [0, s_3 / ft_m], color='grey', linestyle='--')
plt.annotate('Distance Traveled in Climb', xy=(s_trav / nmi_m - 1, 16000), rotation='vertical')
plt.xlabel("Distance Traveled (nmi)")
plt.ylabel("Altitude (ft)")
plt.xlim(0,55)
plt.ylim(0,30200)
plt.show()
