from Constants import *
#Constants
f_con = 5/100                   #-
m_res = 0.30                    #-
a_regression = 0.5617           #- Regression value
b_regression = 1199.7           #- Regression value
total_eff = eta_fuelcell * eta_inverter * eta_wire * eta_EM * eta_prop
m_f_extra = 0.03                #-
fuel_mass_ref = 669             #kg
m_tanks = 1.4 * fuel_mass_ref   #kg
CL_CD = 20
W_P_design = 0.077
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
    mf_mMTO_fraction = 1 - np.exp((-1 * R) / (total_eff * (e_lh2 / g) * CL_CD))
    return mf_mMTO_fraction, R
m_crew_total = m_crew + m_crew_baggage
m_payload = m_pax + m_pax_baggage + m_crew_total
fc_mass = (1 / total_eff * eta_prop * eta_fuelcell) / fc_power_density
em_mass = (1 / eta_EM) / em_power_density
inverter_mass = (1 / eta_EM / eta_wire / eta_inverter) / inverter_power_density
masses_sum = (fc_mass + em_mass + inverter_mass) / W_P_design * g / 10**3

def mtom(oew_ratio):
    # mtom = (m_payload + m_crew_total) / (1 - 2.4 * mf_mMTO("full") * (1 + m_res) - oew_ratio)
    mtom = (m_payload) / (1 - (1) * mf_mMTO("full")[0] * (1 + m_res) - oew_ratio)
    return mtom
def fuel_mass(oew_ratio, range):
    if range == "full":
        m_f = mtom(oew_ratio) * (mf_mMTO(range)[0] * (1 + m_res)) * (1 + m_f_extra)
    elif range == 500:
        m_f = mtom(oew_ratio) * (mf_mMTO(range))[0] * (1 + m_f_extra)
    elif range == 1000:
        m_f = mtom(oew_ratio) * (mf_mMTO(range))[0] * (1 + m_f_extra)
    return m_f

oew_mtom = 12770 / 18287.63
m_mto = mtom(oew_mtom)
m_f = fuel_mass(oew_mtom, "full")
oem = oew_mtom * m_mto #+ m_f * 1.4
m_zf = m_mto - m_f
m_tanks = 1.4 * m_f
print("MTOM:", m_mto)
print("OEM:", oem, oew_mtom)
print("Fuel mass:", fuel_mass(oew_mtom, "full"))
print(mf_mMTO("full"))
print(m_payload)