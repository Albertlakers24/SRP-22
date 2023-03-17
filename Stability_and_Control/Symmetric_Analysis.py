import control
import matplotlib.pyplot as plt
import numpy as np
from Constants.AircraftGeometry import S_w , c_mac_w
from Constants.Stability_Control import CZadot, Cmadot, KY2, CXu, CXa, CZ0, CXq, CZu, CZa, CX0, CZq, Cmu, Cmalpha, Cmq, CXde, CZde, Cmde
from Constants.MissionInputs import rho_0, g, V_cruise
from Constants.Masses_Locations import m_mto

# Changed namings:
S = S_w
c = c_mac_w

# Needed from Other files:
rho0 = rho_0
Lambda = 1
hp0 =1
Temp0 =1
R =1
V0 = V_cruise
alpha0 = 1
theta0 = 0
m= m_mto
# todo continue here


# ---------------------------------------------
# Equations of Symmetric motion  in State-Space
# ---------------------------------------------

def matrix_symmetric(hp0, V0, theta0, m):
    rho = rho0 * pow(((1 + (Lambda * hp0 / Temp0))), (-((g / (Lambda * R)) + 1)))
    W = m * g
    muc = m / (rho * S * c)
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
                     [Cmde] ])

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
    D_s = np.array([ [0],
                     [0],
                     [0],
                     [0] ])

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
    C_1, C_2, C_3 = matrix_symmetric(hp0, V0, theta0, m)
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
        plt.subplot(int(f"{n}{m}1"))  # Make n rows of m graphs, start with plot 1
        plt.plot(time, velocity)
        plt.title("Velocity vs Time")
        plt.ylabel("Velocity [m/s]")
        plt.xlabel("Time [s]")
        plt.grid()

        # Angle of attack vs Time
        # -----------------------
        plt.subplot(int(f"{n}{m}2"))
        plt.plot(time, alpha)
        plt.title("Angle of attack vs Time")
        plt.ylabel("Alpha [rad]")
        plt.xlabel("Time [s]")
        plt.grid()

        # Pitch angle vs Time
        # -------------------
        plt.subplot(int(f"{n}{m}3"))
        plt.plot(time, theta)
        plt.title("Pitch angle vs Time")
        plt.ylabel("Theta [rad]")
        plt.xlabel("Time [s]")
        plt.grid()

        # Pitch rate vs Time
        # ------------------
        plt.subplot(int(f"{n}{m}4"))
        plt.plot(time, pitchrate)
        plt.title("Pitch rate vs Time")
        plt.ylabel("q [rad/s]")
        plt.xlabel("Time [s]")
        plt.grid()

        plt.tight_layout()
        plt.show()

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

        # Step Input Function
        # -------------------
        step_duration = 100                                          # step input duration                       [s/100]
        elevator_deflection = elevator_deflection                    # elevator deflection during step             [rad]
        u_impulse = elevator_deflection * np.ones(step_duration)     # step
        u = np.append(u_impulse, np.zeros(n_steps - step_duration))  # elevator deflection angle [rad] input vs time [s]

        return t, u, x0


plotting = True
# todo determine:
elevator_deflection = -0.005
response_duration = 240
t, u, x0 = step_deflection(elevator_deflection, response_duration)
data = symmetric_forced_response(plotting, t, u, x0, hp0, V0, alpha0, theta0, m)
print(data)