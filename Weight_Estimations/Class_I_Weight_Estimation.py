import numpy as np
from Constants.MissionInputs import h_cruise,V_cruise,R_div,R_norm,m_pax,m_crew,m_pax_baggage,m_crew_baggage,g, E,Aw
from Constants.FlightPerformance_Propulsion import eta_fuelcell,eta_inverter,eta_EM,eta_prop,eta_wire, e_lh2
from Constants.Masses_Locations import W_P_design,W_S_design
from Constants.Aerodynamics import CL_CD_DesCruise
#Constants
f_con = 5/100                   #-
m_res = 0.30                    #-
a_regression = 0.5617           #- Regression value
b_regression = 1199.7           #- Regression value
total_eff = eta_fuelcell * eta_inverter * eta_wire * eta_EM * eta_prop          #[-]
m_f_extra = 0.03                #-
fuel_mass_ref = 669             #kg
m_tanks = fuel_mass_ref/0.2     #kg
CL_CD = CL_CD_DesCruise         #[-]


def mf_mMTO(range):
    if range == "full":
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm + R_lost) * (1 + f_con) + 1.2 * R_div + (E * V_cruise)  # m
    elif range == 500:
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm / 2 + R_lost) * (1 + f_con)  # m
    elif range == 1000:
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm + R_lost) * (1 + f_con)  # m
    mf_mMTO_fraction = 1 - np.exp((-1 * R) / (total_eff * (e_lh2 / g) * CL_CD))  #[-]
    return mf_mMTO_fraction, R

m_crew_total = m_crew + m_crew_baggage                  #[kg]
m_payload = m_pax + m_pax_baggage + m_crew_total        #[kg]

def mtom(oew_ratio):
    mtom = (m_payload + 12402.219695693011) / (1 - mf_mMTO("full")[0] * (1 + m_res))         #[-]
    return mtom
def fuel_mass(oew_ratio, range):
    if range == "full":
        m_f = mtom(oew_ratio) * (mf_mMTO(range)[0] * (1 + m_res)) * (1 + m_f_extra)         #[kg]
    elif range == 500:
        m_f = mtom(oew_ratio) * (mf_mMTO(range))[0] * (1 + m_f_extra)                       #[kg]
    elif range == 1000:
        m_f = mtom(oew_ratio) * (mf_mMTO(range))[0] * (1 + m_f_extra)                       #[kg]
    return m_f

oew_mtom = 12402.219695693011 / 18536.35322811106
m_mto = mtom(oew_mtom)
m_f = fuel_mass(oew_mtom, "full")
oem = oew_mtom * m_mto  #+ m_f * 1.4
m_zf = m_mto - m_f
m_tanks = m_f/0.2
print("MTOM:", m_mto)
print("OEM:", oem, oew_mtom)
print("Fuel mass:", fuel_mass(oew_mtom, "full"))
print("Payload:", m_payload)
print("Tank mass:", m_tanks)
print("ZFM:", m_zf)
# print(mf_mMTO("full"))
P_need = m_mto * g / W_P_design
S = m_mto * g / W_S_design
b = np.sqrt(S*Aw)
print(P_need/1000)
print(S)
print(b)
# #Propeller design
# Dp = 0.55 * (P_need/(1000*4))**(1/4)
# print(Dp)
# span_engine_1 = 3.01/2 + 1 + Dp/2
# print(span_engine_1)
# span_engine_2 = 3.01/2 + 1 + Dp + 0.23 + Dp/2
# print(span_engine_2)
# # span_engine_3 = 3.01/2 + 1 + Dp + 0.23 + Dp + 0.23 + Dp/2
# # print(span_engine_3)
# # print(b/2)