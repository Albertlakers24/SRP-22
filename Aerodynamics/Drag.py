import numpy as np
from Constants.MissionInputs import M_cruise, V_cruise, Calculate_Reynolds_number, V_approach, kv_cruise, kv_TO_Land
from Constants.AircraftGeometry import S_w, c_mac_w, taperw, l_f, l_cp, l_cylindrical, l_tc, bw, delta_flap, c_rw, c_tw
from Constants.Aircraft_Geometry_Drawing import A_base_fuse, Amax_fuse, Amax_nacelle,l_nacelle,upsweep_ang, depth_lg, width_lg, S_lg_front
from Aerodynamics.AeroFunctions import *
from Constants.Empennage_LandingGear import taperh, taperv, xc_m_tail, tc_tail_vt, cmac_v, cmac_h, Sh, Sv, Sweep_LE_HT, Sweep_LE_VT


# PRELIMINARY Values -------------
k = 0.405*10**(-5)                      # K value for production sheet metal
percent_lam_fuse = 0.1                 # % laminar flow on fuselage - FIND!!! based on aircraft type/material
percent_lam_wing = 0.15                  # % laminar flow on wing+tails - FIND!!! based on aircraft type/material
tc_wing = 0.21                          # Average t/c ratio of wing crosssection
xc_m_wing = 0.3483                          # (x/c)_m Position of max thickness - 0.3 for low speed, 0.5 for high speed
sweep_m_wing = Calculate_wingsweep(0, 0.25, xc_m_wing, taperw)  # Sweep at (x/c)_m Wing
sweep_m_Htail = Calculate_wingsweep(np.degrees(Sweep_LE_HT), 0,  xc_m_tail, taperh)  # Sweep at (x/c)_m
sweep_m_Vtail = Calculate_wingsweep(np.degrees(Sweep_LE_VT), 0, xc_m_tail, taperv)  # Sweep at (x/c)_m
d_fuse = np.sqrt(Amax_fuse/np.pi) *2    # Max frontal area of fuselage
l_nc = l_cp
Mach_TO = 0.3                           # Take-off Mach no. - approx
Mach_land = 0.21                        # Landing Mach no. - approx
delta_f_TO = 15                         # Flap deflection - take off
delta_f_LD = delta_flap                  # Flap deflection - landing
CD_LP = 1.05                            # 5-10% due to leakage and Excrescence
V_TO = 55

c_mac = (2/3)*c_rw*((1+taperw+taperw**2)/(1+taperw))
c_range = np.linspace(c_rw, c_tw, 100)
b_range = np.linspace(0,bw/2, 100)
a = np.mean(c_range[np.isclose(b_range, d_fuse/2, rtol=0.1)])
S_unexposed = 2*(c_rw + a)/2*d_fuse/2
# ------------wetted areas +common things!! -----------------
Swet_tailv = 1.05 * 2 * Sv
Swet_tailh= 1.05 * 2 * Sh
Swet_fuse = (np.pi*d_fuse/4) *((1/3/l_nc**2) *((4*l_nc**2 + (d_fuse**2/4))**(1.5) - (d_fuse**3/8)) - d_fuse + 4*l_cylindrical + 2*np.sqrt(l_tc**2 + (d_fuse**2/4))  )
Swet_wing = 1.07* 2 * (S_w - S_unexposed)
Swet_nacelle = 2*np.pi* np.sqrt(Amax_nacelle/np.pi) * l_nacelle #4 * l_nacelle *Amax_nacelle

f_fuse = Calculate_f(l_f, Amax_fuse)
f_nacelle = Calculate_f(l_nacelle, Amax_nacelle)

interference_wing, interference_tail, interference_lg, interference_nacelle, interference_fuse_belly = Calculate_interference_drag(0)

# -------------------------- CRUISE --------------------------

Re_Cruise_wing = Calculate_Reynolds_number(V_cruise, c_mac_w, kv_cruise)
Re_Cruise_Htail= Calculate_Reynolds_number(V_cruise, cmac_h, kv_cruise)
Re_Cruise_VTail =Calculate_Reynolds_number(V_cruise, cmac_v, kv_cruise)
Re_Cruise_Fuselage = Calculate_Reynolds_number(V_cruise, l_f, kv_cruise)
Re_Cruise_nacelle = Calculate_Reynolds_number(V_cruise, l_nacelle, kv_cruise)

Re_DragEst_WingCR = Calculate_Re_DragEst(Re_Cruise_wing,k,c_mac_w)
Re_DragEst_FuselageCR = Calculate_Re_DragEst(Re_Cruise_Fuselage,k,l_f)
Re_DragEst_HTailCR = Calculate_Re_DragEst(Re_Cruise_Htail,k,cmac_h)
Re_DragEst_VTailCR = Calculate_Re_DragEst(Re_Cruise_VTail,k,cmac_v)
Re_DragEst_NacelleCR = Calculate_Re_DragEst(Re_Cruise_nacelle,k,l_nacelle)

Cf_fuselageCR = Calculate_Cf(Re_DragEst_FuselageCR, M_cruise, percent_lam_fuse)
Cf_nacelleCR = Calculate_Cf(Re_DragEst_NacelleCR, M_cruise, percent_lam_fuse)
Cf_WingCR = Calculate_Cf(Re_DragEst_WingCR, M_cruise, percent_lam_wing)
Cf_HTailCR = Calculate_Cf(Re_DragEst_HTailCR, M_cruise, percent_lam_wing)
Cf_VTailCR = Calculate_Cf(Re_DragEst_VTailCR, M_cruise, percent_lam_wing)

FF_wing_CR = Calculate_FF_aerocomp(tc_wing, xc_m_wing, sweep_m_wing, M_cruise)
FF_Htail_CR = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Htail, M_cruise)
FF_Vtail_CR = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Vtail, M_cruise)
FF_fuse = Calculate_FF_fuse(f_fuse)
FF_nacelle = Calculate_FF_nacelle(f_nacelle)

Cd_basedrag_CR = Calculate_basedrag(M_cruise, A_base_fuse)
Cd_upsweep = Calculate_cd_upsweep(upsweep_ang, Amax_fuse)

# -------------------------- TAKE OFF --------------------------
Re_TO_wing = Calculate_Reynolds_number(V_TO, c_mac_w, kv_TO_Land)
Re_TO_Htail= Calculate_Reynolds_number(V_TO, cmac_h, kv_TO_Land)
Re_TO_VTail =Calculate_Reynolds_number(V_TO, cmac_v, kv_TO_Land)
Re_TO_Fuselage = Calculate_Reynolds_number(V_TO, l_f, kv_TO_Land)
Re_TO_nacelle = Calculate_Reynolds_number(V_TO, l_nacelle, kv_cruise)

Re_DragEst_WingTO = Calculate_Re_DragEst(Re_TO_wing,k,c_mac_w)
Re_DragEst_FuselageTO = Calculate_Re_DragEst(Re_TO_Fuselage,k,l_f)
Re_DragEst_HTailTO = Calculate_Re_DragEst(Re_TO_Htail,k,cmac_h)
Re_DragEst_VTailTO = Calculate_Re_DragEst(Re_TO_VTail,k,cmac_v)
Re_DragEst_NacelleTO = Calculate_Re_DragEst(Re_TO_nacelle,k,l_nacelle)

Cf_fuselageTO = Calculate_Cf(Re_DragEst_FuselageTO, Mach_TO, percent_lam_fuse)
Cf_nacelleTO = Calculate_Cf(Re_DragEst_NacelleTO, Mach_TO, percent_lam_fuse)
Cf_WingTO = Calculate_Cf(Re_DragEst_WingTO, Mach_TO, percent_lam_wing)
Cf_HTailTO = Calculate_Cf(Re_DragEst_HTailTO, Mach_TO, percent_lam_wing)
Cf_VTailTO = Calculate_Cf(Re_DragEst_VTailTO, Mach_TO, percent_lam_wing)

FF_wing_TO = Calculate_FF_aerocomp(tc_wing, xc_m_wing, sweep_m_wing, Mach_TO)
FF_Htail_TO = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Htail, Mach_TO)
FF_Vtail_TO = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Vtail, Mach_TO)

Cd_basedrag_TO= Calculate_basedrag(Mach_TO, A_base_fuse)

# -------------------------- LANDING --------------------------

Re_LD_wing = Calculate_Reynolds_number(V_approach, c_mac_w, kv_TO_Land)
Re_LD_Htail= Calculate_Reynolds_number(V_approach, cmac_h, kv_TO_Land)
Re_LD_VTail =Calculate_Reynolds_number(V_approach, cmac_v, kv_TO_Land)
Re_LD_Fuselage = Calculate_Reynolds_number(V_approach, l_f, kv_TO_Land)
Re_LD_nacelle = Calculate_Reynolds_number(V_approach, l_nacelle, kv_cruise)

Re_DragEst_WingLD = Calculate_Re_DragEst(Re_LD_wing,k,c_mac_w)
Re_DragEst_FuselageLD = Calculate_Re_DragEst(Re_LD_Fuselage,k,l_f)
Re_DragEst_HTailLD = Calculate_Re_DragEst(Re_LD_Htail,k,cmac_h)
Re_DragEst_VTailLD = Calculate_Re_DragEst(Re_LD_VTail,k,cmac_v)
Re_DragEst_NacelleLD = Calculate_Re_DragEst(Re_LD_nacelle,k,l_nacelle)

Cf_fuselageLD = Calculate_Cf(Re_DragEst_FuselageLD, Mach_land, percent_lam_fuse)
Cf_nacelleLD = Calculate_Cf(Re_DragEst_NacelleLD, Mach_land, percent_lam_fuse)
Cf_WingLD = Calculate_Cf(Re_DragEst_WingLD, Mach_land, percent_lam_wing)
Cf_HTailLD = Calculate_Cf(Re_DragEst_HTailLD, Mach_land, percent_lam_wing)
Cf_VTailLD = Calculate_Cf(Re_DragEst_VTailLD, Mach_land, percent_lam_wing)

FF_wing_LD = Calculate_FF_aerocomp(tc_wing, xc_m_wing, sweep_m_wing, Mach_land)
FF_Htail_LD = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Htail, Mach_land)
FF_Vtail_LD = Calculate_FF_aerocomp(tc_tail_vt, xc_m_tail, sweep_m_Vtail, Mach_land)

#   -----COMMON -----
Cd_basedrag_LD= Calculate_basedrag(Mach_land, A_base_fuse)

Cd_lg_op = Calculate_cd_lg(depth_lg, width_lg, S_lg_front, 0)
Cd_lg_cl = Calculate_cd_lg(depth_lg, width_lg, S_lg_front, 1)

Cd_flap_TO_slotted= Calculate_cd_flap(1, delta_f_TO)

Cd_flap_LD_slotted= Calculate_cd_flap(1, np.degrees(delta_f_LD))

# -------------------------- CRUISE --------------------------

CD0_wingCR =Cf_WingCR*(FF_wing_CR) *Swet_wing
CD0_tailvCR = (Cf_VTailCR) *(FF_Vtail_CR)*Swet_tailv
CD0_tailhCR = (Cf_HTailCR) *(FF_Htail_CR)*Swet_tailh
CD0_FuselageCR = Cf_fuselageCR* (FF_fuse)*Swet_fuse
CD0_nacelleCR = Cf_nacelleCR *FF_fuse*Swet_nacelle
CD0_miscCR = Cd_upsweep + Cd_basedrag_CR

CD0_CR = CD_LP*((CD0_wingCR*interference_wing + (CD0_tailvCR + CD0_tailhCR)*interference_tail + CD0_FuselageCR*interference_fuse_belly + CD0_nacelleCR*interference_nacelle)/S_w+ CD0_miscCR)

print('CD0 Cruise',CD0_CR)

# -------------------------- TAKE OFF  --------------------------

CD0_wingTO =Cf_WingTO*FF_wing_TO*Swet_wing
CD0_tailvTO = (Cf_VTailTO) *(FF_Vtail_TO)*Swet_tailv
CD0_tailhTO = (Cf_HTailTO) *(FF_Htail_TO)*Swet_tailh
CD0_FuselageTO = Cf_fuselageTO* (FF_fuse)*Swet_fuse
CD0_nacelleTO = Cf_nacelleTO *FF_nacelle*(Swet_nacelle)
CD0_miscTO = Cd_upsweep + Cd_basedrag_TO + Cd_lg_cl + Cd_flap_TO_slotted

CD0_TO = CD_LP*((CD0_wingTO*interference_wing + (CD0_tailvTO + CD0_tailhTO)*interference_tail + CD0_FuselageTO*interference_fuse_belly + CD0_nacelleTO*interference_nacelle)/S_w+ CD0_miscTO)

print('CDO Take off', CD0_TO)

# -------------------------- LANDING  --------------------------

CD0_wingLD =Cf_WingLD*FF_wing_LD*Swet_wing
CD0_tailvLD = (Cf_VTailLD) *(FF_Vtail_LD)*Swet_tailv
CD0_tailhLD = (Cf_HTailLD) *(FF_Htail_LD)*Swet_tailh
CD0_FuselageLD = Cf_fuselageLD* (FF_fuse)*Swet_fuse
CD0_nacelleLD = Cf_nacelleLD*FF_fuse*Swet_nacelle
CD0_miscLD = Cd_upsweep + Cd_basedrag_LD + Cd_lg_cl + Cd_flap_LD_slotted

CD0_LD = CD_LP*((CD0_wingLD*interference_wing + (CD0_tailvLD + CD0_tailhLD)*interference_tail + CD0_FuselageLD*interference_fuse_belly + CD0_nacelleLD*interference_nacelle)/S_w+ CD0_miscLD)

print('CD0 Landing', CD0_LD)

