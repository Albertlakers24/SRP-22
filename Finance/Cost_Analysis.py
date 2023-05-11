import numpy as np
from matplotlib import pyplot as plt
from Constants.MissionInputs import V_cruise, kts_m_s, lbs_kg,kW_hp
from Constants.Masses_Locations import m_mto,m_oem
def CEF(year):
    inflation_rate = 6.31752 + 0.104415*(year-2017)
    return inflation_rate
V_max = V_cruise * (1/kts_m_s)
W_to = m_mto * (1/lbs_kg)
W_E = m_oem * (1/lbs_kg)
W_eng = 447
W_fc = 80 * 11
W_A = (W_E - W_eng - W_fc) * (1/lbs_kg)
P_max_hp = 2500 * kW_hp
P_max = 2500
Cef_1989 = 3.02
Cef_1990 = 3.05
Cef_2022 = CEF(2022)
Cef_2035 = CEF(2035)
Cef_2014 = CEF(2014)
###Research, Development, Test and Evaluation cost ###
def Cost_RDTE(C_aedr,C_dstr,C_ftar,C_ftor,C_tsfr,C_pror,C_finr):
    c_RDTE = (C_aedr + C_dstr + C_ftar + C_ftor)/ (1- (C_tsfr + C_pror + C_finr))
    return np.round(c_RDTE,0)
def W_ampr(mass):
    W_ampr = 10 ** (0.1936 + 0.8645 * (np.log10(mass)))
    return np.round(W_ampr,0)
def C_aedr(W_ampr,n_rdte, f_diff,f_cad,r_er):
    Mhr_aedr = 0.0396 * W_ampr**0.791 * V_max**1.526 * n_rdte**0.183 * f_diff * f_cad
    C_aedr = Mhr_aedr * r_er
    return np.round(C_aedr,0)
def C_dstr(W_ampr,CEF,f_diff,n_rdte):
    C_dstr = 0.008325*W_ampr**0.873 * V_max**1.890 * n_rdte**0.346 * CEF * f_diff
    return np.round(C_dstr,0)
def C_ftar(C_ear,C_manr,C_matr,C_toolr,C_qcr):
    C_ftar = C_ear + C_manr + C_matr + C_toolr + C_qcr
    return C_ftar
def C_ear(c_eng, c_prop,c_avionics,n_rdte,n_st,c_tank):
    C_ear = (c_eng + c_prop + c_avionics + c_tank) * (n_rdte - n_st)
    return np.round(C_ear)
def C_prop(power,n_p,cef):
    C_pr = 10 ** (0.7746 + 1.1432*np.log10(power)) * n_p * Cef_2022/ cef
    return np.round(C_pr,0)
def C_eng(c_em,c_fc,p_max,cfe_2022,cef,type):
    if type == 1:
        C_eng = (c_em*p_max*cef/cfe_2022) + (c_fc*p_max*cef/cfe_2022)
    if type ==2 :
        C_eng = (c_em*p_max*cef/cfe_2022) + ((c_fc*0.5)*p_max*cef/cfe_2022)
    if type ==3 :
        C_eng = (c_em*p_max*cef/cfe_2022) + (c_fc* p_max)
    return np.round(C_eng,0)
def C_avionics(mass,c_avionics_perc,cef):
    AMP_est = 10 ** (1.1846 + 1.2625*np.log10(mass))
    C_avionics = AMP_est * c_avionics_perc * cef/Cef_1989
    return np.round(C_avionics,0)
def C_tank(S_tank,rho_cb,c_cb,thickness,c_cb_m):
    C_tank = (S_tank * thickness * rho_cb) * (c_cb + c_cb_m)
    return C_tank
def C_manr(W_ampr,n_rdte,f_diff,Rmr):
    Mhr_manr = 28.984 * W_ampr**0.740 * V_max**0.543 * n_rdte**0.524 * f_diff
    C_manr = Mhr_manr * Rmr
    return np.round(C_manr,0)
def C_matr(f_mat,W_ampr,n_rdte,CEF):
    C_matr = 37.632 * f_mat * W_ampr**0.689 * V_max**0.624 * n_rdte**0.792 * CEF
    return np.round(C_matr,0)
def C_toolr(W_ampr,n_rdte,n_rr,f_diff,Rtr):
    Mhr_tool = 4.0127*W_ampr**0.764 * V_max**0.899 * n_rdte**0.178 * n_rr**0.066 * f_diff
    C_toolr = Mhr_tool * Rtr
    return np.round(C_toolr,0)
def C_qcr(C_manr):
    C_qcr = 0.13 * C_manr
    return np.round(C_qcr,0)
def C_ftor(W_ampr,n_rdte,n_st,f_diff,f_obs,cef):
    C_ftor = 0.001244*W_ampr**1.160 * V_max**1.371 * (n_rdte - n_st)**1.281*cef*f_diff*f_obs
    return np.round(C_ftor,0)
Rer_1989 = 62.0
Rer_2035 = Rer_1989 * Cef_2035/Cef_1989
Rmr_1989 = 34.44
Rmr_2035 = Rmr_1989 * Cef_2035/Cef_1989
Rtr_1989 = 43.06
Rtr_2035 = Rtr_1989 * Cef_2035/Cef_1989
N_rdte = 6
N_st = 3
N_rr = 0.33
F_diff = 1.7
F_cad = 1.12000
F_mat = 1.0
F_obs = 1.0
N_p = 12
P_per_prop_hp = P_max_hp / N_p
C_tsfr = 0.2
C_pror = 0.10
C_finr = 0.10
C_avionics_perc = 0.10
C_em = 240 * 1.5
C_fc = 5000
C_fc_opimistic = 640
C_CF = 90 * Cef_2035/Cef_2022
C_CF_manufacturing = 220000/1000 * Cef_2035/ Cef_2014
S_tank = 23.5
tank_t = 5.6 / 1000
rho_CF = 2000
C_aedr = C_aedr(W_ampr(W_to),N_rdte,F_diff,F_cad,Rer_2035)
C_dstr = C_dstr(W_ampr(W_to),Cef_2035,F_diff,N_rdte)
C_eng1 = C_eng(C_em,C_fc,P_max,Cef_2022,Cef_2035,1)
C_eng2 = C_eng(C_em,C_fc,P_max,Cef_2022,Cef_2035,2)
C_eng3 = C_eng(C_em,C_fc_opimistic,P_max,Cef_2022,Cef_2035,3)
C_avionics = C_avionics(W_to,C_avionics_perc,Cef_2035)
C_prop = C_prop(P_per_prop_hp,N_p,Cef_2035)
C_tank = C_tank(S_tank, rho_CF, C_CF, tank_t, C_CF_manufacturing)
C_ear1 = C_ear(C_eng1,C_prop,C_avionics,N_rdte,N_st,C_tank)
C_ear2 = C_ear(C_eng2,C_prop,C_avionics,N_rdte,N_st,C_tank)
C_ear3 = C_ear(C_eng3,C_prop,C_avionics,N_rdte,N_st,C_tank)
C_manr = C_manr(W_ampr(W_to),N_rdte,F_diff,Rmr_2035)
C_matr = C_matr(F_mat,W_ampr(W_to),N_rdte,Cef_2035)
C_toolr = C_toolr(W_ampr(W_to),N_rdte,N_rr,F_diff,Rtr_2035)
C_qcr = C_qcr(C_manr)
C_ftar1 = C_ftar(C_ear1,C_manr,C_matr,C_toolr,C_qcr)
C_ftar2 = C_ftar(C_ear2,C_manr,C_matr,C_toolr,C_qcr)
C_ftar3 = C_ftar(C_ear3,C_manr,C_matr,C_toolr,C_qcr)
C_ftor = C_ftor(W_ampr(W_to),N_rdte,N_st,F_diff,F_obs,Cef_2035)
C_RDTE1 = Cost_RDTE(C_aedr,C_dstr,C_ftar1,C_ftor,C_tsfr,C_pror,C_finr)
C_RDTE2 = Cost_RDTE(C_aedr,C_dstr,C_ftar2,C_ftor,C_tsfr,C_pror,C_finr)
C_RDTE3 = Cost_RDTE(C_aedr,C_dstr,C_ftar3,C_ftor,C_tsfr,C_pror,C_finr)
# labels = ["Airframe engineering and design","Development support and testing", "Flight test airplane", "Flight test operations","Profit","Finance","Test, simulation facilities"]
# plt.pie([C_aedr/C_RDTE1,C_dstr/C_RDTE1,C_ftar1/C_RDTE1,C_ftor/C_RDTE1,C_pror,C_finr,C_tsfr],autopct='%1.1f%%',pctdistance=1.2,textprops={'fontsize': 12})
# plt.title("Research, Development and Testing cost for no fuel cell price reduction")
# plt.legend(labels,bbox_to_anchor=(1, 1), loc=2, frameon=False,fontsize=12)
# plt.show()
# plt.pie([C_aedr/C_RDTE2,C_dstr/C_RDTE2,C_ftar2/C_RDTE2,C_ftor/C_RDTE2,C_pror,C_finr,C_tsfr],autopct='%1.1f%%',pctdistance=1.2,textprops={'fontsize': 12})
# plt.title("Research, Development and Testing cost for 50% fuel cell price reduction")
# plt.legend(labels,bbox_to_anchor=(1, 1), loc=2, frameon=False,fontsize=12)
# plt.show()
# all_masses = np.array([C_aedr,C_dstr,C_ftar3,C_ftor,C_pror * C_RDTE3,C_finr * C_RDTE3,C_tsfr * C_RDTE3]) / 10 **6
# labels = [f'{l}, {s:1.1f} Million USD' for l, s in zip(labels, all_masses)]
# plt.pie([C_aedr/C_RDTE3,C_dstr/C_RDTE3,C_ftar3/C_RDTE3,C_ftor/C_RDTE3,C_pror,C_finr,C_tsfr])
#plt.title("Research, Development and Testing cost for optimistic fuel cell price")
# plt.legend(labels,bbox_to_anchor=(0.95, 0.75),fontsize=12)
# plt.show()
# print(C_RDTE1/10**6)
# print(C_RDTE2/10**6)
# print(C_RDTE3/10**6)

###Manufacturing and Acquisition cost ###
AEP_list1 = []
Manufacturing_cost_lst1 = []
AEP_list2 = []
Manufacturing_cost_lst2 = []
AEP_list3 = []
Manufacturing_cost_lst3 = []
C_aedm_lst = []
C_eam_lst1 = []
C_eam_lst2 = []
C_eam_lst3 = []
C_int_lst = []
C_manm_lst = []
C_matm_lst = []
C_toolm_lst = []
C_qcm_lst = []
C_tankm_lst = []
C_ftom_lst = []
for N_m in range(1,2000,1):
    N_program = N_m + N_rdte
    F_int = 1000
    N_pax = 48
    R_em = Rer_2035
    R_mm = Rmr_2035
    R_tm = Rtr_2035
    N_rm = 3
    T_pft = 10
    F_ftoh = 4.0
    C_finm = 0.10
    F_prom = 0.2
    C_opshr = 2310 * Cef_2035/Cef_2022
    def C_ACQ(c_man, c_pro):
        C_ACQ = c_man + c_pro
        return np.round(C_ACQ,-6)
    def AEP(c_acq,c_rdte,n_m):
        AEP = (c_acq + c_rdte)/n_m
        return np.round(AEP,0)
    def C_man(c_aedm, c_apcm, c_ftom, c_finm):
        C_man = (c_aedm + c_apcm + c_ftom) / (1- c_finm)
        return np.round(C_man,0)
    def C_aedm(mhr_aed_program,r_em,c_aedr):
        C_aedm = mhr_aed_program * r_em - c_aedr
        return np.round(C_aedm,0)
    def Mhr_aed_program(W_ampr, n_program,f_diff,f_cad):
        Mhr_aed_program = 0.0396*W_ampr**0.791 * V_max**1.526 * n_program**0.183 * f_diff * f_cad
        return Mhr_aed_program
    def C_apcm(c_eam, c_intm, c_manm, c_matm, c_toolm, c_qcm, c_tankm):
        C_apcm = c_eam + c_intm + c_manm + c_matm + c_toolm + c_qcm + c_tankm
        return np.round(C_apcm,0)
    def C_eam(c_eng,c_prop,c_avionics,n_m):
        C_eam = (c_eng + c_prop +c_avionics) * n_m
        return np.round(C_eam,0)
    def C_intm(f_int,n_pax,n_m,cef_1990,cef):
        C_intm = f_int * n_pax * n_m * cef/cef_1990
        return np.round(C_intm,0)
    def C_manm(w_ampr,n_program,f_diff,r_mm,c_manr):
        mhr_man_program = 28.984*w_ampr**0.740 * V_max**0.543 * n_program**0.524 * f_diff
        C_manm = mhr_man_program * r_mm - c_manr
        return np.round(C_manm,0)
    def C_matm(f_mat,w_ampr,n_program,cef,c_matr):
        c_mat_program = 37.632*f_mat * w_ampr**0.689 * V_max**0.624 * n_program**0.792 * cef
        C_matm = c_mat_program - c_matr
        return np.round(C_matm,0)
    def C_toolm(w_ampr, n_program,n_rm,f_diff,r_tm,c_toolr):
        mhr_tool_program = 4.0127*w_ampr**0.764 * V_max**0.899 * n_program**0.178 * n_rm**0.066 * f_diff
        C_toolm = mhr_tool_program * r_tm - c_toolr
        return np.round(C_toolm,0)
    def C_qcm(c_manm):
        C_qcm = 0.13 * c_manm
        return np.round(C_qcm,0)
    def C_tankm(c_tank,n_m):
        C_tankm = c_tank * n_m
        return C_tankm
    def C_ftom(n_m,c_opshr, t_pft, f_ftoh):
        C_ftom = n_m * c_opshr * t_pft * f_ftoh
        return np.round(C_ftom,0)
    def C_pro(f_prom, c_man):
        C_pro = f_prom * c_man
        return np.round(C_pro,0)
    C_aedm = C_aedm(Mhr_aed_program(W_ampr(W_to),N_program,F_diff,F_cad),R_em,C_aedr)
    C_eam1 = C_eam(C_eng1,C_prop,C_avionics,N_m)
    C_eam2 = C_eam(C_eng2, C_prop, C_avionics, N_m)
    C_eam3 = C_eam(C_eng3, C_prop, C_avionics, N_m)
    C_intm = C_intm(F_int,N_pax,N_m,Cef_1990,Cef_2035)
    C_manm = C_manm(W_ampr(W_to),N_program,F_diff,R_mm,C_manr)
    C_matm = C_matm(F_mat,W_ampr(W_to),N_program,Cef_2035,C_matr)
    C_toolm = C_toolm(W_ampr(W_to),N_program,N_rm,F_diff,R_tm,C_toolr)
    C_qcm = C_qcm(C_manm)
    C_tankm = C_tankm(C_tank,N_m)
    C_apcm1 = C_apcm(C_eam1,C_intm,C_manm,C_matm,C_toolm,C_qcm,C_tankm)
    C_apcm2 = C_apcm(C_eam2, C_intm, C_manm, C_matm, C_toolm, C_qcm,C_tank)
    C_apcm3 = C_apcm(C_eam3, C_intm, C_manm, C_matm, C_toolm, C_qcm,C_tank)
    C_ftom = C_ftom(N_m,C_opshr,T_pft,F_ftoh)
    C_man1 = C_man(C_aedm,C_apcm1,C_ftom,C_finm)
    C_man2 = C_man(C_aedm, C_apcm2, C_ftom, C_finm)
    C_man3 = C_man(C_aedm, C_apcm3, C_ftom, C_finm)
    C_pro1 = C_pro(F_prom, C_man1)
    C_pro2 = C_pro(F_prom, C_man2)
    C_pro3 = C_pro(F_prom, C_man3)
    C_ACQ1 = C_ACQ(C_man1,C_pro1)
    C_ACQ2 = C_ACQ(C_man2, C_pro2)
    C_ACQ3 = C_ACQ(C_man3, C_pro3)
    AEP1 = AEP(C_ACQ1,C_RDTE1,N_m)
    AEP2 = AEP(C_ACQ2, C_RDTE2, N_m)
    AEP3 = AEP(C_ACQ3, C_RDTE3,N_m)
    C_aedm_lst.append(C_aedm/10**6)
    C_ftom_lst.append(C_ftom/10**6)
    C_eam_lst1.append(C_eam1/10**6)
    C_eam_lst2.append(C_eam2/10**6)
    C_eam_lst3.append(C_eam3/10**6)
    C_int_lst.append(C_intm/10**6)
    C_manm_lst.append(C_manm/10**6)
    C_matm_lst.append(C_matm/10**6)
    C_toolm_lst.append(C_toolm/10**6)
    C_qcm_lst.append(C_qcm/10**6)
    C_tankm_lst.append(C_tankm/10**6)
    AEP_list1.append(AEP1/10**6)
    Manufacturing_cost_lst1.append(C_man1/ 10**6)
    AEP_list2.append(AEP2/10 ** 6)
    Manufacturing_cost_lst2.append(C_man2 / 10**6)
    AEP_list3.append(AEP3/10 ** 6)
    Manufacturing_cost_lst3.append(C_man3 / 10**6)
# plt.plot(np.arange(1,2000,1),AEP_list1,label = "Without Fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),AEP_list2, label = "With 50% Fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),AEP_list3, label = "With optimistic Fuel cell price")
# plt.ylim(0,100)
# plt.xlabel("Number of aircraft sold",fontsize=12)
# plt.ylabel("USD in million",fontsize=12)
# plt.legend(fontsize=12)
# plt.show()
# plt.plot(np.arange(1,2000,1),AEP_list1[999]*np.arange(1,2000,1), label = "Revenue without fuel cell price reduction")
# plt.hlines(C_RDTE1/10**6,0,2000,label = "RDTE cost without fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),Manufacturing_cost_lst1 + C_RDTE1/10**6, label = "Total cost without fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),AEP_list2[999]*np.arange(1,2000,1), label = "Revenue with 50% fuel cell price reduction")
# plt.hlines(C_RDTE2/10**6,0,2000,label = "RDTE cost with 50% fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),Manufacturing_cost_lst2 + C_RDTE2/10**6, label = "Total cost with 50% fuel cell price reduction")
# plt.plot(np.arange(1,2000,1),AEP_list3[999]*np.arange(1,2000,1), label = "Revenue optimistic fuel cell price")
# plt.hlines(C_RDTE3/10**6,0,2000,label = "RDTE cost with optimistic fuel cell price")
# plt.plot(np.arange(1,2000,1),Manufacturing_cost_lst3 + C_RDTE3/10**6, label = "Total cost with optimistic fuel cell price")
# plt.xlabel("Number of aircraft sold")
# plt.ylabel("USD in million")
# plt.xlim(0,1000)
# plt.ylim(0,50000)
# plt.legend()
# plt.show()
#
# labels = ["Airframe engineering and design","Engine and avionics", "Interior","Manufacturing labor","Manufacturing material","Tooling","Tank","Quality control","Production flight test operations","Finance"]
# plt.pie([C_aedm_lst[999]/Manufacturing_cost_lst1[999],C_eam_lst1[999]/Manufacturing_cost_lst1[999],C_int_lst[999]/Manufacturing_cost_lst1[999],C_manm_lst[999]/Manufacturing_cost_lst1[999],
#          C_matm_lst[999]/Manufacturing_cost_lst1[999],C_toolm_lst[999]/Manufacturing_cost_lst1[999],C_tankm_lst[999]/Manufacturing_cost_lst1[999],C_qcm_lst[999]/Manufacturing_cost_lst1[999],
#          C_ftom_lst[999]/Manufacturing_cost_lst1[999],C_finm],autopct='%1.1f%%',pctdistance=1.1,textprops={'rotation_mode':'anchor','fontsize': 12})
# plt.title("Manufacturing cost for without fuel cell price reduction")
# plt.legend(labels,bbox_to_anchor=(1, 1), loc=2, frameon=False,fontsize=12)
# plt.show()
# plt.pie([C_aedm_lst[999]/Manufacturing_cost_lst2[999],C_eam_lst2[999]/Manufacturing_cost_lst2[999],C_int_lst[999]/Manufacturing_cost_lst2[999],C_manm_lst[999]/Manufacturing_cost_lst2[999],
#          C_matm_lst[999]/Manufacturing_cost_lst2[999],C_toolm_lst[999]/Manufacturing_cost_lst2[999],C_tankm_lst[999]/Manufacturing_cost_lst2[999],C_qcm_lst[999]/Manufacturing_cost_lst2[999],
#          C_ftom_lst[999]/Manufacturing_cost_lst2[999],C_finm],autopct='%1.1f%%',pctdistance=1.1,radius = 1)
# plt.title("Manufacturing cost for 50% fuel cell price reduction")
# plt.legend(labels,bbox_to_anchor=(1, 1), loc=2, frameon=False)
# plt.show()
# plt.pie([C_aedm_lst[999]/Manufacturing_cost_lst3[999],C_eam_lst3[999]/Manufacturing_cost_lst3[999],C_int_lst[999]/Manufacturing_cost_lst3[999],C_manm_lst[999]/Manufacturing_cost_lst3[999],
#          C_matm_lst[999]/Manufacturing_cost_lst3[999],C_toolm_lst[999]/Manufacturing_cost_lst3[999],C_tankm_lst[999]/Manufacturing_cost_lst3[999],C_qcm_lst[999]/Manufacturing_cost_lst3[999],
#          C_ftom_lst[999]/Manufacturing_cost_lst3[999],C_finm],autopct='%1.1f%%',pctdistance=1.1,radius = 1)
# plt.title("Manufacturing cost for optimistic fuel cell price")
# plt.legend(labels,bbox_to_anchor=(1, 1), loc=2, frameon=False)
# plt.show()

# #Break even point calculation
# for i in range(0,2000):
#     Revenue = i * AEP_list1[999]
#     if Revenue >= (C_RDTE1/10**6 + Manufacturing_cost_lst1[i]):
#         print("It takes",i,"of aircraft sold to break even")
#         ROI = (AEP_list1[999]*1000 - (C_RDTE1/10**6 + Manufacturing_cost_lst1[999])) / (C_RDTE1/10**6 + Manufacturing_cost_lst1[999]) * 100
#         print("The ROI is:",np.round(ROI,2))
#         print("The AEP is:", AEP_list1[999])
#         print("The AEP in 2022:", AEP_list1[999]*Cef_2022/Cef_2035)
#         print("The Manufacturing cost per aircraft is:", Manufacturing_cost_lst1[999]/1000)
#         break
#     else:
#         continue
#
# for i in range(0,2000):
#     Revenue = i * AEP_list2[999]
#     if Revenue >= (C_RDTE2/10**6 + Manufacturing_cost_lst2[i]):
#         print("It takes",i,"of aircraft sold to break even")
#         ROI = (AEP_list2[999]*1000 - (C_RDTE2 / 10 ** 6 + Manufacturing_cost_lst2[999])) / (C_RDTE2 / 10 ** 6 + Manufacturing_cost_lst2[999]) *100
#         print("The ROI is:", np.round(ROI,2))
#         print("The AEP is:",AEP_list2[999])
#         print("The AEP in 2022:", AEP_list2[999]*Cef_2022/Cef_2035)
#         print("The Manufacturing cost per aircraft is:", Manufacturing_cost_lst2[999]/1000)
#         break
#     else:
#         continue
#
# for i in range(0,2000):
#     Revenue = i * AEP_list3[999]
#     if Revenue >= (C_RDTE3/10**6 + Manufacturing_cost_lst3[i]):
#         print("It takes",i,"of aircraft sold to break even")
#         ROI = (AEP_list3[999]*1000 - (C_RDTE3 / 10 ** 6 + Manufacturing_cost_lst3[999])) / (C_RDTE3 / 10 ** 6 + Manufacturing_cost_lst3[999]) * 100
#         print("The ROI is:",np.round(ROI,2))
#         print("The AEP is:", AEP_list3[999])
#         print("The AEP in 2022:", AEP_list3[999]*Cef_2022/Cef_2035)
#         print("The Manufacturing cost per aircraft is:", Manufacturing_cost_lst3[999]/1000)
#         break
#     else:
#         continue
# print(Manufacturing_cost_lst1[999])
# print(Manufacturing_cost_lst2[999])
# print(Manufacturing_cost_lst3[999])

### Operating Cost
R_bl1 = 500
R_bl2 = 1000
t_cl = 20.9/60
R_cl = 66.5
R_de = 110
t_de = 0.5462
AHj = 850
T_efj = 11.0 * Cef_2035/Cef_1990
Sal_j_cap = 52000 * Cef_2035/Cef_1989 #20000 * Cef_2035/Cef_1989
Sal_j_fo = 21000 * Cef_2035/Cef_1989 #11000 * Cef_2035/Cef_1989
W_f_500nmi = 241
W_f_full = 683
F_ins_hull = (0.005+0.030)/2
R_l_ap = 10.0 * Cef_2035/Cef_1989 #(10.0+22.60)/2 * Cef_2035/Cef_1989  #10.0 * Cef_2035/Cef_1989
R_l_eng = 10.0 * Cef_2035/Cef_1989  #(10.0+22.60)/2 * Cef_2035/Cef_1989  # 10.0 * Cef_2035/Cef_1989
F_p = 2                     #Price of hydrogen $/kg could potentially go down to 1$/kg
f_amb_lab = 1.2
f_amb_mat = 0.55
F_dap = 0.85
DP_ap = 10
F_deng = 0.85
DP_eng = 10
F_dprp = 0.85
DP_prp = 7
F_dav = 1
DP_av = 5
F_dapsp = 0.85
DP_apsp = 10
F_dengsp = 0.85
DP_engsp = 7
F_dfc = 0.85
DP_fc = 5
C_rt = 0.001 + 10**(-8) * W_to
DOC_fin = 0.07
C_ins = 0.02
H_em = 5000
EP = P_max * C_em * Cef_2035/Cef_2022
AFP1 = AEP_list1[999] - C_eng1
AFP2 = AEP_list2[999] - C_eng2
AFP3 = AEP_list3[999] - C_eng3
def R_bl_annual(v_bl,u_annual_bl):
    R_bl_annual = v_bl * u_annual_bl
    return R_bl_annual
def V_bl(r_bl,t_bl):
    V_bl = r_bl/ t_bl
    return V_bl
def T_bl(t_gm,t_cl,t_cr,t_de):
    T_bl = t_gm + t_cl + t_cr + t_de
    return T_bl
def T_gm(mass):
    T_gm = 0.51 * 10**(-6) * mass + 0.125
    return T_gm
def T_cr(r_bl,r_cl,r_de,r_man,v_max):
    T_cr = (1.01*r_bl - r_cl - r_de - r_man)/v_max
    return T_cr
def T_man(mass):
    T_man = 0.25 * 10 **(-6) * mass + 0.025
    return T_man
def R_man(t_man,type):
    if type == 1:
        R_man = 250 * t_man
    if type == 2:
        R_man = V_cruise/kts_m_s * t_man
    return R_man
def U_ann_bl(t_bl):
    U_ann_bl = 10**3 * (3.4546*t_bl + 2.994 - (12.289*t_bl**2 - 5.56626*t_bl +8.964)**(1/2))
    return U_ann_bl
def V_flt(v_cr,t_cr,t_cl,t_de):
    V_flt = v_cr * (t_cr / (t_cl + t_cr + t_de))
    return V_flt
def t_flt(t_cl,t_cr,t_de):
    t_flt = t_cl + t_cr + t_de
    return t_flt
def DOC(doc_flt,doc_maint,doc_depr,doc_lnr,doc_fin,c_rt,c_ins):
    DOC = (doc_flt + doc_maint + doc_depr + doc_lnr) / (1 - doc_fin - c_rt - c_ins)
    return DOC
def DOC_flt(c_crew,c_pol):
    DOC_flt = c_crew + c_pol
    return DOC_flt
def C_crew(n_cap,n_cp,v_bl,ahj,salj_cap,salj_fo,tefj):
    C_cap = n_cap * (1.26/v_bl) * (salj_cap/ahj) + tefj/v_bl
    C_cp = n_cp * (1.26/v_bl) * (salj_fo/ahj) + tefj/v_bl
    return C_cap + C_cp
def C_pol(w_f,r_bl,f_p):
    C_pol = w_f/r_bl * f_p * 1.05
    return C_pol
def DOC_maint(c_lab_ap,c_lab_eng,c_mat_ap,c_mat_eng,c_amb):
    DOC_maint = c_lab_ap + c_lab_eng + c_mat_ap + c_mat_eng + c_amb
    return DOC_maint
def C_lab_ap(mass,r_l_ap,v_bl):
    C_lab_ap = 1.03 * (3.0 + 0.067*mass/1000) * r_l_ap / v_bl
    return C_lab_ap
def C_lab_eng(p_max,n_eng,r_l_eng,v_bl):
    C_lab_eng = 1.03 * 1.3 * n_eng *((0.4956 + 0.0532*(p_max/n_eng/1000))*(1100/H_em)+0.1) * r_l_eng / v_bl * 0.5 # the 0.5 can change if found a reduction with better electrical maintainance
    return C_lab_eng
def C_mat_ap(CEF,v_bl,AEP,N_e,EP):
    c_mat_apblhr = 30 * (CEF/Cef_1989) + 0.79*10**(-5)*(AEP-(N_e*EP))
    C_mat_ap = 1.03 * c_mat_apblhr / v_bl
    return C_mat_ap
def C_mat_eng(N_e,v_bl,EP):
    K_hem = 0.021*(H_em/100) + 0.769
    C_mat_eng_blhr = ((5.43*10**(-5) * EP)*1.5 - 0.47) * (1/K_hem) * 0.5 # the 0.5 can change if found a reduction with better electrical maintainance
    C_mat_eng = 1.03 *1.3 * N_e * C_mat_eng_blhr / v_bl
    return C_mat_eng
def C_amb(v_bl, mass,p_max,n_eng,AEP,EP,CEF):
    Mhr_map_bl = 3 + 0.067*mass/1000
    Mhr_eng_bl = ((0.4956 + 0.0532*(p_max/n_eng/1000))*(1100/H_em)+0.1)
    c_mat_apblhr = 30 * (CEF/Cef_1989) + 0.79*10**(-5)*(AEP-(n_eng*EP))
    K_hem = 0.021*(H_em/100) + 0.769
    C_mat_eng_blhr = ((5.43*10**(-5) * EP)*1.5 - 0.47) * (1/K_hem) * 0.5
    C_amb = 1.03* ((f_amb_lab * (Mhr_map_bl*R_l_ap + n_eng*Mhr_eng_bl*R_l_eng) ) + (f_amb_mat*(c_mat_apblhr + n_eng*C_mat_eng_blhr))) / v_bl
    return C_amb
def DOC_depr(C_dap, C_deng, C_dprp, C_dva, C_dapsp, C_dengsp,C_dfc):
    DOC_depr = C_dap + C_deng + C_dprp + C_dva + C_dapsp + C_dengsp + C_dfc
    return DOC_depr
def C_dap(AEP,N_e,EP,N_p,PP,ASP,FCP,U_annual_bl,V_bl):
    C_dap = (F_dap*(AEP - N_e*EP - N_p*PP - ASP - FCP))/(DP_ap * U_annual_bl * V_bl)
    return C_dap
def C_deng(n_e,EP,U_annual_bl,V_bl):
    C_deng = F_deng * n_e * EP / (DP_eng * U_annual_bl*V_bl)
    return C_deng
def C_dprp(n_p,PP,U_annual_bl,V_bl):
    C_dprp = F_dprp * n_p * PP / (DP_prp * U_annual_bl*V_bl)
    return C_dprp
def C_dav(ASP,U_annual_bl,V_bl):
    C_dav = F_dav * ASP / (DP_av * U_annual_bl*V_bl)
    return C_dav
def C_dapsp(AEP,n_e,EP,FCP,U_annual_bl,V_bl):
    C_dapsp = F_dapsp*0.1 *(AEP - n_e * EP - FCP) / (DP_apsp * U_annual_bl * V_bl)
    return C_dapsp
def C_dengsp(n_e,EP,U_annual_bl,V_bl):
    C_dengsp = F_dengsp*0.5*(n_e*EP*1.5)/(DP_engsp*U_annual_bl*V_bl)
    return C_dengsp
def C_dfc(FCP,U_annual_bl,V_bl):
    C_dfc = F_dfc * FCP / (DP_fc * U_annual_bl*V_bl)
    return C_dfc
def DOC_lnr(C_lf,C_nf):
    DOC_lnr = C_lf + C_nf
    return DOC_lnr
def C_lf(mass,V_bl,t_bl):
    C_lf = 0.002*mass / (V_bl*t_bl)
    return C_lf
def C_nf(V_bl,t_bl):
    C_nf = 10*Cef_2035/Cef_1989 /(V_bl*t_bl)
    return C_nf
def IOC (DOC):
    IOC = 0.52/0.48 * DOC
    return IOC

T_gm = T_gm(W_to)
T_man = T_man(W_to)
R_man = R_man(T_man,1)
t_cr = T_cr(R_bl1,R_cl,R_de,R_man,V_max)
t_flt = t_flt(t_cl,t_cr,t_de)
t_bl = T_bl(T_gm,t_cl, t_cr, t_de)
V_bl = V_bl(R_bl1,t_bl)
U_annual_bl = U_ann_bl(t_bl)
R_bl_annual = R_bl_annual(V_bl,U_annual_bl)
C_crew = C_crew(1,1,V_bl,AHj,Sal_j_cap,Sal_j_fo,T_efj)
C_pol = C_pol(W_f_500nmi,R_bl1,F_p)
DOC_flt = DOC_flt(C_crew,C_pol)
C_lab_ap = C_lab_ap(W_A,R_l_ap,V_bl)
C_lab_eng = C_lab_eng(P_max_hp,4,R_l_eng,V_bl)
C_mat_ap1 = C_mat_ap(Cef_2035,V_bl,AEP_list1[999],1,C_eng1)
C_mat_ap2 = C_mat_ap(Cef_2035,V_bl,AEP_list2[999],1,C_eng2)
C_mat_ap3 = C_mat_ap(Cef_2035,V_bl,AEP_list3[999],1,C_eng3)
C_mat_eng = C_mat_eng(1,V_bl,EP)
C_amb1 = C_amb(V_bl,W_A,P_max_hp,1,AEP_list1[999],C_eng1,Cef_2035)
C_amb2 = C_amb(V_bl,W_A,P_max_hp,1,AEP_list2[999],C_eng2,Cef_2035)
C_amb3 = C_amb(V_bl,W_A,P_max_hp,1,AEP_list3[999],C_eng3,Cef_2035)
DOC_maint1 = DOC_maint(C_lab_ap,C_lab_eng,C_mat_ap1,C_mat_eng,C_amb1)
DOC_maint2 = DOC_maint(C_lab_ap,C_lab_eng,C_mat_ap2,C_mat_eng,C_amb2)
DOC_maint3 = DOC_maint(C_lab_ap,C_lab_eng,C_mat_ap3,C_mat_eng,C_amb3)
C_dap1 = C_dap(AEP_list1[999],1,EP,1,C_prop,C_avionics,P_max*C_fc*Cef_2035/Cef_2022,U_annual_bl, V_bl)
C_dap2 = C_dap(AEP_list2[999],1,EP,1,C_prop,C_avionics,P_max*C_fc*Cef_2035/Cef_2022*0.5,U_annual_bl, V_bl)
C_dap3 = C_dap(AEP_list3[999],1,EP,1,C_prop,C_avionics,P_max*C_fc_opimistic,U_annual_bl, V_bl)
C_deng = C_deng(1,EP,U_annual_bl, V_bl)
C_dprp = C_dprp(1,C_prop,U_annual_bl, V_bl)
C_dav = C_dav(C_avionics, U_annual_bl, V_bl)
C_dapsp1 = C_dapsp(AEP_list1[999],1,EP,P_max*C_fc*Cef_2035/Cef_2022,U_annual_bl, V_bl)
C_dapsp2 = C_dapsp(AEP_list2[999],1,EP,P_max*C_fc*Cef_2035/Cef_2022*0.5 ,U_annual_bl, V_bl)
C_dapsp3 = C_dapsp(AEP_list3[999],1,EP,P_max*C_fc_opimistic,U_annual_bl, V_bl)
C_dengsp = C_dengsp(1,EP,U_annual_bl, V_bl)
C_dfc1 = C_dfc(P_max*C_fc*Cef_2035/Cef_2022,U_annual_bl, V_bl)
C_dfc2 = C_dfc(P_max*C_fc*Cef_2035/Cef_2022*0.5,U_annual_bl, V_bl)
C_dfc3 = C_dfc(P_max*C_fc_opimistic,U_annual_bl, V_bl)
DOC_depr1 = DOC_depr(C_dap1,C_deng,C_dprp,C_dav,C_dapsp1,C_dengsp, C_dfc1)
DOC_depr2 = DOC_depr(C_dap2,C_deng,C_dprp,C_dav,C_dapsp2,C_dengsp, C_dfc2)
DOC_depr3 = DOC_depr(C_dap3,C_deng,C_dprp,C_dav,C_dapsp3,C_dengsp, C_dfc3)
C_lf = C_lf(W_to,V_bl,t_bl)
C_nf = C_nf(V_bl, t_bl)
DOC_lnr = DOC_lnr(C_lf,C_nf)
DOC1 = DOC(DOC_flt,DOC_maint1,DOC_depr1,DOC_lnr,DOC_fin,C_rt,C_ins)
DOC2 = DOC(DOC_flt,DOC_maint2,DOC_depr2,DOC_lnr,DOC_fin,C_rt,C_ins)
DOC3 = DOC(DOC_flt,DOC_maint3,DOC_depr3,DOC_lnr,DOC_fin,C_rt,C_ins)
# all_cost = np.array([DOC_flt,DOC_maint3,DOC_depr3,DOC_lnr,DOC_fin * DOC3,C_rt* DOC3,C_ins* DOC3])*V_bl
# all_cost_label = ["Flying Cost", "Maintenance Cost", "Depreciation Cost", "Landing Fess", "Finance Cost","Tax Cost"]
# labels = [f'{l}, {s:1.1f} USD/hr' for l, s in zip(all_cost_label, all_cost)]
# plt.pie(all_cost)
# plt.title("Direct Operating Cost for optimistic fuel cell price")
# plt.legend(labels,bbox_to_anchor=(0.95, 0.75),fontsize=12)
# plt.show()

Total_OC1 = DOC1 + IOC(DOC1)
Total_OC2 = DOC2 + IOC(DOC2)
Total_OC3 = DOC3 + IOC(DOC3)
C_ops_airline1 = Total_OC1 * R_bl_annual * 25 * 10
C_ops_airline2 = Total_OC2 * R_bl_annual * 25 * 10
C_ops_airline3 = Total_OC3 * R_bl_annual * 25 * 10
PAX = 48
#print(C_crew * V_bl *Cef_2022/Cef_2035)
# print(DOC1/PAX)
# print(DOC2/PAX)
# print(DOC3/PAX)
# print(DOC1*Cef_2022/Cef_2035/PAX)
# print(DOC2*Cef_2022/Cef_2035/PAX)
# print(DOC3*Cef_2022/Cef_2035/PAX)
# print(Total_OC1/PAX)
# print(Total_OC2/PAX)
# print(Total_OC3/PAX)
# print(C_ops_airline1/10**6)
# print(C_ops_airline2/10**6)
# print(C_ops_airline3/10**6)