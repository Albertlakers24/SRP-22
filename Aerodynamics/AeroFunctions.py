from Constants.MissionInputs import M_cruise
from Constants.AircraftGeometry import S_w, taper, c_mac, t_c_ratio, Aw
import pandas as pd
import numpy as np

CLmax_Clmax = 0.9          # Wing/Airfoil CL -- From graph - ADSEE II, Lecture 2, p.19
eta_airfoil = 0.95          # Airfoil efficiency -- ADSEE Slides

def Calculate_beta(mach):
    beta = np.sqrt(1 - (mach ** 2))
    return beta

def Calculate_wingsweep(sweep1, position1, position2, taper):
    sweep = np.degrees(np.arctan(np.tan(np.radians(sweep1)) - (4/Aw*(position2 - position1) *(1-taper)/(1+taper))))
    return sweep

def Calculate_CL_alpha(beta, sweep):
    CL_alpha = (2*np.pi*Aw)/ (2 + np.sqrt(4 + ((Aw *beta/eta_airfoil)**2 * (1 + ((np.tan(np.radians(sweep)))**2)/(beta**2)))))
    return CL_alpha

def Calculate_CL_max(Cl_max, delta_Clmax):
    CL_max = CLmax_Clmax*Cl_max + delta_Clmax
    return CL_max

def Calculate_alpha_stall(CL_max, CL_alpha, alpha_0L, delta_alpha_CLMax):
    alpha_stall = (CL_max/CL_alpha) #+ alpha_0L + delta_alpha_CLMax
    return alpha_stall

def Calculate_alpha_trim(CL_des, CL_alpha, Alpha0):
    alpha_trim = CL_des/CL_alpha + Alpha0
    return alpha_trim

def Making_labels(filename):
    with open(filename) as f:
        first_line = f.readline().split()
    return first_line

def Reading_data(filename, labels):
    data_dictionary = pd.read_csv(filename,
                               sep="\s+",
                               skiprows=1,
                               usecols=np.arange(0, len(labels), 1),
                               names=labels)
    return data_dictionary


def Calculate_Cf(Re, mach, percent_lam):
    Cf_laminar = 1.328/np.sqrt(Re)
    Cf_turb = 0.455/((np.log10(Re))**(2.58) * ((1 + (0.144 * mach**2)))**(0.65)  )
    Cf_tot = (1-percent_lam)*Cf_turb + percent_lam*Cf_laminar
    return Cf_tot


def Calculate_Re_DragEst(Re, k, l):
    cutoff_Re= 38.21*(l/k)**(1.053)              # Check that (l/k) is corrected substituted!!
    Re_DragEst = min(Re, cutoff_Re)
    return Re_DragEst

def Calculate_FF_aerocomp(t_c_avg, x_c_m, sweep_m, mach):     # WING, TAIL, STRUT, PYLON
    FF = ( 1 + (0.6/x_c_m*t_c_avg) + (100*t_c_avg**4)) * (1.34 * mach**0.18* (np.cos(np.radians(sweep_m)))**0.28)
    return FF
    # (x/c)_m ----- position of max thickness
    # (t_c) ----> Average t/c ratio
    # sweep_m = sweep at this position

def Calculate_f(l, Amax):
    f = l/np.sqrt((4/np.pi*Amax))
    return f
    # l --- length of part
    # Amax ---- max frontal area

def Calculate_FF_fuse(f):       # FUSELAGE + SMOOTH CANOPY
    FF = 1 + (60/f**3) + (f/400)
    return FF

def Calculate_FF_nacelle(f):
    FF = 1 +(0.35/f)
    return FF

def Calculate_cd_upsweep(upsweep_ang, Amax):
    cd_upsweep = 3.83 * (np.radians(upsweep_ang))**(2.5) * Amax/S_w                    # CROSS CHECK S DIVISION HERE
    return cd_upsweep

def Calculate_basedrag(mach, A_base):       # Abase ---- fuselage base area
    cd_basedrag = (0.319 + (0.419*(mach - 0.161)**2)) * A_base/ S_w      # CROSS CHECK S DIVISION HERE
    return cd_basedrag

def Calculate_cd_lg(d_tire, w_tire, S_A, wheel_well):
    if wheel_well == 0:
        delta_cds = 0.05328 * np.exp(5.615*S_A/d_tire/w_tire)
    elif wheel_well == 1:
        delta_cds = 0.04955 * np.exp(5.615 * S_A/d_tire/w_tire)
    # else:
    #     print("Error, Invalid wheel well type")
    cd_lg = delta_cds * d_tire*w_tire/S_w
    return cd_lg

def Calculate_cd_flap(flap_type, c_f, S_flap, delta_f):
    if flap_type == 0:
        F_flap = 0.0144
    if flap_type == 1:
        F_flap = 0.0074
    # else:
    #     print("Error, Invalid flap type")

    cd_flap = F_flap * c_f/c_mac * S_flap/S_w *(delta_f - 10)
    return cd_flap