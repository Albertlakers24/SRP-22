import control
import numpy as np
import matplotlib.pyplot as plt
from Constants.AircraftGeometry import S_w , c_mac_w, bw
from Constants.Stability_Control import muc,Cnp,Clr,Clbeta,Cnbeta,Cybeta, Cn_r, KZ2,KY2, Izz, Clp, mub, KX2, CZadot, Cmadot, CXu, CXa, CZ0, CXq, CZu, CZa, CX0, CZq, Cmu, Cmalpha, Cmq, CXde, CZde, Cmdelta
from Constants.MissionInputs import rho_0, g, V_cruise, ISA_calculator, h_cruise,dt_cruise
from Constants.Masses_Locations import m_mto
from Constants.Aerodynamics import CL_DesCruise

# Needed from Other files:
V0 = V_cruise
rho = ISA_calculator(h_cruise,dt_cruise)[2]

# ---------------------------------------------
# Eigenvalues Symmetric Motion
# ---------------------------------------------

# Short Period
# ------------
Cnbeta=0.058125750001392446
A = 2*muc * KY2 * (2*muc - CZadot)
B = -2*muc * KY2 * CZa - (2*muc + CZq) * Cmadot - (2*muc - CZadot) * Cmq
C = CZa * Cmq - (2*muc + CZq) * Cmalpha
eigenvalue1_SP = (complex(-B / (2*A), np.sqrt(4*A*C - B**2) / (2*A))) * (V0/c_mac_w)
eigenvalue2_SP = (complex(-B / (2*A), -np.sqrt(4*A*C - B**2) / (2*A))) * (V0/c_mac_w)

Real_sp = eigenvalue1_SP.real
Imag_sp = eigenvalue1_SP.imag
DampingRatio_sp = - (Real_sp/np.sqrt(Real_sp**2+Imag_sp**2))

# Phugoid Motion
# --------------
A = 2*muc * (CZa * Cmq - 2*muc * Cmalpha)
B = 2*muc * (CXu * Cmalpha - Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa)
C = CZ0 * (Cmu * CZa - CZu * Cmalpha)
eigenvalue1_FM = (complex(- B / (2*A), np.sqrt(4*A*C - B**2) / (2*A))) * (V0/c_mac_w)
eigenvalue2_FM = (complex(- B / (2*A), -np.sqrt(4*A*C - B**2) / (2*A))) * (V0/c_mac_w)

Real_ph = eigenvalue2_FM.real
Imag_ph = eigenvalue2_FM.imag
DampingRatio_ph = - (Real_ph/np.sqrt(Real_ph**2+Imag_ph**2))

# Aperiodic Rolling Motion
eigenvalue_ARM = (Clp/(4*mub*KX2))*(V0/bw)

Period_r_ARM = -1/eigenvalue_ARM

# Dutch Roll
a = 8*mub ** 2 * KZ2
b1 = -2*mub * (Cn_r + 2 * KZ2 * Cybeta)
c = 4*mub * Cnbeta + Cybeta * Cn_r

eigenvalue1_DR = (complex(-b1/ (2*a), np.sqrt(4*a*c - b1**2) / (2*a))) * (V0/bw)
eigenvalue2_DR = (complex(-b1 / (2*a), -np.sqrt(4*a*c - b1**2) / (2*a))) * (V0/bw)

Real_eigenvalue =eigenvalue1_DR.real
Imag_eigenvalue =eigenvalue1_DR.imag
DampingRatio_DR = -Real_eigenvalue/np.sqrt(Real_eigenvalue**2+Imag_eigenvalue**2)
ohm_nd_DR = np.sqrt(Cnbeta*0.2*rho*V_cruise**2*S_w*bw/Izz)
Damping_ohmnd_DR =DampingRatio_DR*ohm_nd_DR

# Aperiodic Spiral Motion
eigenvalue_ASM= (2*CL_DesCruise*(Clbeta*Cn_r-Cnbeta*Clr))/(Clp*(Cybeta*Cn_r+4*mub*Cnbeta)-Cnp*(Cybeta*Clr+4*mub*Clbeta))

Period2s_ASM = -np.log(2)/eigenvalue_ASM*bw/V_cruise

# Add Labels and sizing nicely for Symmetric Cases
plt.xlabel("Real")
plt.ylabel("Imaginary")

"""Graphs for dynamic stability"""
# plt.plot(eigenvalue1_SP.real, eigenvalue1_SP.imag, marker="x", color="r", linewidth=1)        # short period
# plt.annotate("  \u03BB\u2085", (eigenvalue1_SP.real, eigenvalue1_SP.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue2_SP.real, eigenvalue2_SP.imag, marker="x", color="r", linewidth=1)        # short period
# plt.annotate("  \u03BB\u2086", (eigenvalue2_SP.real, eigenvalue2_SP.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue1_FM.real, eigenvalue1_FM.imag, marker="x", color="r", linewidth=1)        # phugoid
# plt.annotate("  \u03BB\u2087", (eigenvalue1_FM.real, eigenvalue1_FM.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue2_FM.real, eigenvalue2_FM.imag, marker="x", color="r", linewidth=1)        # phugoid
# plt.annotate("  \u03BB\u2084", (eigenvalue2_FM.real, eigenvalue2_FM.imag), ha='left', va='top', fontsize=20)
# plt.xlim(-4,0.5)
# plt.ylim(-5,5)
# plt.grid(True)
# plt.axhline(0, color='k', linewidth=0.5)
# plt.axvline(0, color='k', linewidth=0.5)
# plt.show()
#
# # Add Labels and sizing nicely for Asymmetric Cases
# plt.xlabel("Real")
# plt.ylabel("Imaginary")
# plt.plot(eigenvalue_ARM.real, eigenvalue_ARM.imag, marker="x", color="r", linewidth=1)        # Aperiodic Rolling Motion
# plt.annotate("  \u03BB\u2081", (eigenvalue_ARM.real, eigenvalue_ARM.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue1_DR.real, eigenvalue1_DR.imag, marker="x", color="r", linewidth=1)        # Dutch Roll
# plt.annotate("  \u03BB\u2082", (eigenvalue1_DR.real, eigenvalue1_DR.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue2_DR.real, eigenvalue2_DR.imag, marker="x", color="r", linewidth=1)        # Dutch Roll
# plt.annotate("  \u03BB\u2083", (eigenvalue2_DR.real, eigenvalue2_DR.imag), ha='left', fontsize=20)
# plt.plot(eigenvalue_ASM.real, eigenvalue_ASM.imag, marker="x", color="r", linewidth=1)        # Aperiodic Spiral Motion
# plt.annotate("  \u03BB\u2088", (eigenvalue_ASM.real, eigenvalue_ASM.imag), ha='left', fontsize=20)
# plt.xlim(-4,1)
# plt.ylim(-0.6,0.6)
# plt.grid(True)
# plt.axhline(0, color='k', linewidth=0.5)
# plt.axvline(0, color='k', linewidth=0.5)
# plt.show()

print("Symmetric")
print("---------------Short Period----------------")
if DampingRatio_sp >0.3:
    if DampingRatio_sp <2.00:
        print("Damping ratio =", DampingRatio_sp)
    else:
        print("Damping ratio requirement is NOT met:")
        print("Damping ratio =", DampingRatio_sp)

print("---------------Phugoid Motion----------------")
if DampingRatio_ph > 0.0:
    print("Damping ratio =", DampingRatio_ph)
else:
    print("Damping ratio requirement is NOT met:")
    print("Damping ratio =", DampingRatio_ph)

print("Asymmetric")
print("---------------Aperiodic Roll----------------")
if Period_r_ARM>1.4:
    print("Period_r requirement is not met:")
print("Period_r =", Period_r_ARM)

print("---------------Dutch Roll----------------")
if DampingRatio_DR<0.08:
    print("Dampring ratio requirement is not met:")
print("Damping ratio =", DampingRatio_DR)
if Damping_ohmnd_DR<0.15:
    print("Dampring ratio*ohm_nd requirement is not met:")
print("Damping ratio*ohm_nd =", Damping_ohmnd_DR)
if ohm_nd_DR<0.4:
    print("ohm_nd requirement is NOT met:")
print("ohm_nd =", ohm_nd_DR)

print("---------------Spiral Motion----------------")
if Period2s_ASM>8:
    print("Period2s =", Period2s_ASM)
else:
    print("Period2s requirement is not met:")
    print("Period2s =", Period2s_ASM)

# ---------------------------------------------
# Equations of Symmetric motion  in State-Space
# ---------------------------------------------
"""
alpha0 = 1
theta0 = 0

def matrix_symmetric(hp0, V0, alpha0, theta0, m):
    W = m_mto * g
    CX0 = W * np.sin(theta0) / (0.5 * rho * V0 ** 2 * S)
    CZ0 = -W * np.cos(theta0) / (0.5 * rho * V0 ** 2 * S)

    C_1 = np.array([ [-2*muc * (c/V0), 0, 0, 0],
                     [0, CZadot -2*muc * (c/V0), 0, 0],
                     [0, 0, -c/V0, 0],
                     [0, Cmadot * (c/V0), 0, -2*muc * KY2 * (c/V0)] ])

    C_2 = np.array([ [CXu, CXa, CZ0, CXq],
                     [CZu, CZa, -CX0, CZq + 2*muc],
                     [0, 0, 0, 1],
                     [Cmu, Cmalpha, 0, Cmq] ])

    C_3 = np.array([ [CXde],
                     [CZde],
                     [0],
                     [Cmdelta]])
    return C_1, C_2, C_3

def matrix_symmetric_system(C_1, C_2, C_3):
    # State Matrix
    # ------------
    A_s = -np.matmul(np.linalg.inv(C_1), C_2)
    # Input Matrix
    # ------------
    B_s = -np.matmul(np.linalg.inv(C_1), C_3)
    # Output Matrix
    # ------------
    C_s = np.array([ [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1] ])
    # Feedthrough / Feedforward Matrix
    # --------------------------------
    D_s = np.array([[0],
                     [0],
                     [0],
                     [0]])

    return A_s, B_s, C_s, D_s

def find_eigenvalues(A, print_eigenvalues):
    eigenvalues = np.linalg.eigvals(A)
    eigenvalue1_SP = eigenvalues[0]
    eigenvalue2_SP = eigenvalues[1]
    eigenvalue1_FM = eigenvalues[2]
    eigenvalue2_FM = eigenvalues[3]

    if print_eigenvalues == True:
        print(f"                                           \n"
              f"Numerical Model Symmetrical Eigenvalues    \n"
              f"----------------------------------------   \n"
              f"Short Period: {eigenvalue1_SP}             \n"
              f"Short Period: {eigenvalue2_SP}             \n"
              f"Phugoid Motion: {eigenvalue1_FM}           \n"
              f"Phugoid Motion: {eigenvalue2_FM}             ")

    return eigenvalue1_SP, eigenvalue2_SP, eigenvalue1_FM, eigenvalue2_FM

def symmetric_forced_response(plotting, t, u, x0 ,hp0, V0, alpha0, theta0, m):

    # Symmetrical Motion State-Space System
    # -------------------------------------
    C_1, C_2, C_3 = matrix_symmetric(hp0, V0, alpha0, theta0, m)
    A, B, C, D = matrix_symmetric_system(C_1, C_2, C_3)
    ssmodel_s = control.ss(A, B, C, D)

    # Eigenvalues for space-state A matrix
    # ------------------------------------
    if plotting == True:
        print_eigenvalues = True
    else:
        print_eigenvalues = False
    eigenvalue1_SP, eigenvalue2_SP, eigenvalue1_FM, eigenvalue2_FM = find_eigenvalues(A, print_eigenvalues)

    # Symmetrical state-space model response
    # --------------------------------------
    time, y = control.forced_response(ssmodel_s, t, u, x0)
    print("t&y=", time,y)
    # substitutions from dimentionless coefficients
    # ---------------------------------------------
    velocity = V0 + V0 * y[0]            # dimentionless velocity to True Airspeed     [m/s]
    pitchrate = (V0/c) * y[3]            # dimentionless pitch rate to pitch rate    [rad/s]

    # substitutions from stability axis sytem to body axis system
    # -----------------------------------------------------------
    alpha = alpha0 + y[1]          # angle of attack   [rad]
    theta = theta0 + y[2]          # pitch angle       [rad]

    # data formating
    # --------------
    data = dict()

    data['outputs'] = dict()
    data['outputs']['velocity'] = velocity
    data['outputs']['alpha'] = alpha
    data['outputs']['theta'] = theta
    data['outputs']['pitchrate'] = pitchrate
    data['outputs']['time'] = time

    data['eigenvalues'] = dict()
    data['eigenvalues']['short_period'] = eigenvalue1_SP, eigenvalue2_SP
    data['eigenvalues']['phugoid_motion'] = eigenvalue1_FM, eigenvalue2_FM

    # Response Plots for Symmetrical Motion
    # -------------------------------------
    if plotting == True:
        n = 2
        m = 2
        # Velocity vs Time
        # ----------------
        # plt.subplot(int(f"{n}{m}1"))  # Make n rows of m graphs, start with plot 1
        # plt.plot(time, velocity)
        # plt.title("Velocity vs Time")
        # plt.ylabel("Velocity [m/s]")
        # plt.xlabel("Time [s]")
        # #plt.xlim(0,15)                          # For short motion
        # plt.grid()
        #
        # # Angle of attack vs Time
        # # -----------------------
        # plt.subplot(int(f"{n}{m}2"))
        # plt.plot(time, alpha)
        # plt.title("Angle of attack vs Time")
        # plt.ylabel("Alpha [rad]")
        # plt.xlabel("Time [s]")
        # #plt.xlim(0, 15)                         # For short motion
        # plt.grid()
        #
        # # Pitch angle vs Time
        # # -------------------
        # plt.subplot(int(f"{n}{m}3"))
        # plt.plot(time, theta)
        # plt.title("Pitch angle vs Time")
        # plt.ylabel("Theta [rad]")
        # plt.xlabel("Time [s]")
        # #plt.xlim(0, 15)                         # For short motion
        # plt.grid()
        #
        # # Pitch rate vs Time
        # # ------------------
        # plt.subplot(int(f"{n}{m}4"))
        # plt.plot(time, pitchrate)
        # plt.title("Pitch rate vs Time")
        # plt.ylabel("q [rad/s]")
        # plt.xlabel("Time [s]")
        # #plt.xlim(0, 15)                         # For Short motion
        # plt.grid()
        #
        # plt.tight_layout()
        # #plt.show()

    return data

def step_deflection(elevator_deflection, response_duration):
        # Time array
        # ----------
        start = 0                          # begining of time plot     [s]
        stop = response_duration           # end of time plot          [s]
        step = 0.01                        # step size for time plot   [s]

        t = np.arange(start, stop + step, step)
        n_steps = int((stop - start) / step + 1)          # total number of time steps

        # State-space initial state of input vector
        # -----------------------------------------
        x0 = np.array([[0],           # initial velocity dimentionless parameter         [-]
                       [0],           # initial angle of attack                        [rad]
                       [0],           # initial pitch angle                            [rad]
                       [0]])          # initial pitch rate dimentionless coefficient     [-]

        # Step Input Function for Phugoid
        # -------------------
        step_duration = 100                                          # step input duration                       [s/100]
        elevator_deflection = elevator_deflection                    # elevator deflection during step             [rad]
        u_impulse = elevator_deflection * np.ones(step_duration)     # step
        u = np.append(u_impulse, np.zeros(n_steps - step_duration))  # elevator deflection angle [rad] input vs time [s]
        # u = [u,u]                                                    # For trimtab
        # Step Input Function for Short Period
        # step_duration_short = 1500                                     # step input duration                        [s/100]
        # start_short = 0                                                # beginning of time plot                     [s]
        # stop_short = response_duration                                 # end of time plot                           [s]
        # step_short = 0.01                                              # step size for time plot                    [s]
        # t = np.arange(start_short, stop_short + step_short, step_short)
        # n_steps_short = int((stop_short - start_short) / step_short + 1)                     # total number of time steps
        #
        # u_impulse_short = elevator_deflection_short*np.ones(step_duration_short)
        # u = np.append(u_impulse_short, np.zeros(n_steps_short-step_duration_short))
        return t, u, x0

plotting = True
elevator_deflection = -0.005                                           # Small change in elevator deflection for Phugoid  [rad]
elevator_deflection_short = -0.005                                   # Change in elevator deflection for short period   [rad]
response_duration = 240                                                # Originally at 240 sec                            [sec]

t, u, x0 = step_deflection(elevator_deflection, response_duration)
data = symmetric_forced_response(plotting, t, u, x0, hp0, V0, alpha0, theta0, m)
print(data)

# Symmetric Motion Eigenvalues
# ----------------------------
print(f"                                            \n"
      f"Analytical Model Symmetrical Eigenvalues    \n"
      f"----------------------------------------    \n"
      f"Short Period: {eigenvalue1_SP}              \n"
      f"Short Period: {eigenvalue2_SP}              \n"
      f"Phugoid Motion: {eigenvalue1_FM}            \n"
      f"Phugoid Motion: {eigenvalue2_FM}              ")

print("Damping Ratio Short Period=", DampingRatio_sp)
print("Damping Ratio Phugoid=", DampingRatio_ph)

# For Assymetric Motion
plt.xlabel("Real")
plt.ylabel("Imaginary")
plt.plot(-1.35, 0, marker="x", color="r", linewidth=1)             # Aperiodic Roll
plt.annotate("  \u03BB\u2085", (-1.35, 0), ha='left')
plt.plot(-0.201, 1.226, marker="x", color="r", linewidth=1)         # Dutch Roll
plt.annotate("  \u03BB\u2086", (-0.201, 1.226), ha='left')
plt.plot(-0.201, -1.226, marker="x", color="r", linewidth=1)        # Dutch Roll
plt.annotate("  \u03BB\u2087", (-0.201, -1.226), ha='left')
plt.plot(0.00364, 0, marker="x", color="r", linewidth=1)             # Spiral Motion
plt.annotate("  \u03BB\u2088", (0.00364, 0), ha='left')
plt.xlim(-1.8,0.3)
plt.ylim(-1.5,1.5)
plt.grid(True)
plt.axhline(0, color='k', linewidth=0.5)
plt.axvline(0, color='k', linewidth=0.5)
plt.show()
"""