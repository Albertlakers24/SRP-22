import numpy as np
from Cit_par22 import *
import control.matlab as ml
import matplotlib.pyplot as plt

# makes matrices C1, C2 and C3 as described in the V&V plan

#Oscialltion Period
def oscillation_period(eta):
    P = (2 * np.pi) / eta
    return P
# Half amplitude Period
def half_amplitude(epsilon):
    T_1_2 = np.log(1/2)/ epsilon
    return T_1_2
#Damping
def damping(epsilon,eta):
    damp = -epsilon/ np.sqrt(epsilon**2 + eta**2)
    return damp
#Undamped frequency
def undamped_frequency(epsilon, eta):
    omega = np.sqrt(eta**2 + epsilon**2)
    return omega

def matrices_C(motion_type: str, cp: object):
    '''
    Converting the Equations of motions to matrices C1, C2 and C3.
    So that the matrices could then be used later to transform into steady state format.
    The matrices transformation could be found in the V&V plan.
    '''
    if motion_type == "s":
        c = cp.c
        V0 = cp.V0
        muc = cp.muc
        C1 = (c/V0) * np.array([[-2 * muc, 0, 0, 0],
                        [0, (cp.CZadot - 2 * muc), 0, 0],
                        [0, 0, -1, 0],
                        [0, cp.Cmadot, 0, -2 * muc * cp.KY2]])  # do we indeed use V0 for V?

        C2 = np.array([[cp.CXu, cp.CXa, cp.CZ0, cp.CXq],
                        [cp.CZu, cp.CZa, -cp.CX0, (cp.CZq + 2 * muc)],
                        [0, 0, 0, 1],
                        [cp.Cmu, cp.Cma, 0, cp.Cmq]])

        C3 = np.array([[cp.CXde],
                        [cp.CZde],
                        [0],
                        [cp.Cmde]])
    elif motion_type == "a":
        mub = cp.mub
        b = cp.b
        V0 = cp.V0
        C1 = (b/V0) *np.array([[cp.CYbdot - 2*cp.mub, 0, 0, 0],
                               [0, -1/2, 0, 0],
                               [0, 0, (-4* mub*cp.KX2), (4* mub*cp.KXZ)],
                               [cp.Cnbdot, 0, (4* mub*cp.KXZ), (-4* mub*cp.KZ2)]])  # do we indeed use V0 for V?
        C2 = np.array([[cp.CYb, cp.CL, cp.CYp, (cp.CYr - 4*mub)],
                       [0, 0, 1, 0],
                       [cp.Clb, 0, cp.Clp, cp.Clr],
                       [cp.Cnb, 0, cp.Cnp, cp.Cnr]])
        C3 = np.array([[cp.CYda, cp.CYdr],
                        [0, 0],
                        [cp.Clda, cp.Cldr],
                        [cp.Cnda, cp.Cndr]])
    return C1, C2, C3


# makes steady state matrices as described in the V&V plan
def ss_matrices(C1: np.ndarray, C2: np.ndarray, C3: np.ndarray, motion_type: str):
    '''
    Matrices C1,C2 and C3 into steady state format and returning steady state matrices A,B,C,D.
    The matrices transformation could be found in the V&V plan.
    '''
    A = np.matmul(-1 * np.linalg.inv(C1), C2)
    B = np.matmul(-1 * np.linalg.inv(C1), C3)
    C = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    if motion_type == "s":
        D = np.array([[0],
                      [0],
                      [0],
                      [0]])
    elif motion_type == "a":
        D = np.array([[0, 0],
                      [0, 0 ],
                      [0, 0],
                      [0, 0]])
    return A, B, C, D


# Symmetric motion
def sym_state_space(cp):
    EOM_matrices_sym = matrices_C("s",cp)
    ss_matrices_sym = ss_matrices(EOM_matrices_sym[0], EOM_matrices_sym[1], EOM_matrices_sym[2], "s")
    sym_motion = ml.ss(ss_matrices_sym[0], ss_matrices_sym[1], ss_matrices_sym[2], ss_matrices_sym[3])
    # poles,zeroes = ml.pzmap(sym_motion)
    # plt.annotate('$\lambda_3$',(poles[0].real + 0.1,poles[0].imag))
    # plt.annotate('$\lambda_4$', (poles[1].real + 0.1, poles[1].imag))
    # plt.annotate('$\lambda_1$', (poles[2].real + 0.1, poles[2].imag))
    # plt.annotate('$\lambda_2$', (poles[3].real+ 0.1, poles[3].imag))
    # plt.xlim(-3,0.5)
    # plt.show()
    # Phuogid
    # print("Phugoid Oscillation Period:",oscillation_period(np.abs(poles[2].imag)))
    # print("Phugoid Half Amplitude:",half_amplitude(poles[2].real))
    # print("Phugoid Damp:",damping(poles[2].real,np.abs(poles[2].imag)))
    # print("Phugoid undamped frequency:",undamped_frequency(poles[2].real,np.abs(poles[2].imag)))
    #
    # # Short Period
    # print("Short Period Oscillation Period:", oscillation_period(np.abs(poles[0].imag)))
    # print("Short Period Half Amplitude:", half_amplitude(poles[0].real))
    # print("Short Period Damp:", damping(poles[0].real,np.abs(poles[0].imag)))
    # print("Short Period undamped frequency:", undamped_frequency(poles[0].real,np.abs(poles[0].imag)))
    return sym_motion

# Asymmetric motion
def asym_state_space(cp):
    EOM_matrices_asym = matrices_C("a",cp)
    ss_matrices_asym = ss_matrices(EOM_matrices_asym[0],EOM_matrices_asym[1],EOM_matrices_asym[2],"a")
    asym_motion = ml.ss(ss_matrices_asym[0],ss_matrices_asym[1],ss_matrices_asym[2],ss_matrices_asym[3])
    poles, zeroes = ml.pzmap(asym_motion)
    # plt.annotate('$\lambda_1$', (poles[0].real + 0.1, poles[0].imag))
    # plt.annotate('$\lambda_3$', (poles[1].real + 0.1, poles[1].imag))
    # plt.annotate('$\lambda_4$', (poles[2].real + 0.1, poles[2].imag))
    # plt.annotate('$\lambda_2$', (poles[3].real + 0.1, poles[3].imag))
    # plt.xlim(-3, 0.5)
    # plt.show()

    # Aperiodic Roll
    # print("Aperiodic Roll Amplitude:", half_amplitude(poles[0].real))
    # print("Aperiodic Roll Damp:", damping(poles[0].real, poles[0].imag))
    # print("Aperiodic Roll undamped frequency:", undamped_frequency(poles[0].real, np.abs(poles[0].imag)))
    #
    # # Dutch Roll
    # print("Dutch Roll Oscillation Period:", oscillation_period(np.abs(poles[2].imag)))
    # print("Dutch Roll Half Amplitude:", half_amplitude(poles[2].real))
    # print("Dutch Roll Damp:", damping(poles[2].real, np.abs(poles[2].imag)))
    # print("Dutch Roll undamped frequency:", undamped_frequency(poles[2].real, np.abs(poles[2].imag)))
    # #
    # # Spiral
    # print("Spiral Half Amplitude:", half_amplitude(poles[3].real))
    # print("Spiral Damp:", damping(poles[3].real, np.abs(poles[3].imag)))
    # print("Spiral undamped frequency:", undamped_frequency(poles[3].real, np.abs(poles[3].imag)))
    return asym_motion


def x_vector_sym(u, alpha, theta, q):
    return np.array([[u],
                     [alpha],
                     [theta],
                     [q]])


def x_vector_asym(beta, phi, p, r):
    return np.array([[beta],
                     [phi],
                     [p],
                     [r]])


def u_vector_sym(delta_e):
    return delta_e


def u_vector_asym(delta_a, delta_r):
    return np.array([delta_a,
                     delta_r])



# t = np.arange(0, 20.01, 0.1)
# # Phugoid control input
# d_e = np.array([0., 0.1, 0.3, 0.5, 0.4, 0.2, 0.0])
# d_e = np.append(d_e, np.zeros(194))
# d_e = np.transpose(d_e)
# u_phugoid = u_vector_sym(d_e)
# x0 = x_vector_sym(u=150, alpha=2.5, theta=0.5, q=11250)
# # print(ss_matrices_sym)
# #y_out, t_out, x_out = ml.lsim(sym_motion, u_phugoid, t, x0)
# y_out,t = ml.impulse(sym_motion)
#print('xout',x_out)
# print('yout', y_out)
# aoa = []
# for i in range(len(y_out)):
#     aoa.append(y_out[i][1])
# plt.plot(t, aoa)
# plt.show()
# print()
# Short period control input

# A-periodic roll control input

# Dutch roll control input
# def broodje():
#     t_spi = np.arange(0, 120.1, 0.1)
#     X0 =
#
# def csv_data_into_np(filename):
#     return np.genfromtxt(filename)


# Phugoid
def phugoid(start_index, end_index, cp):
    t = np.arange(0, 150.1, 0.1)
    U_real = np.radians(np.genfromtxt('FlightData/delta_e.csv')[start_index:end_index])
    V0 = np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index]
    theta0 = np.radians(np.genfromtxt('FlightData/Ahrs1_Pitch.csv')[start_index])
    alpha0 = np.radians(np.genfromtxt('FlightData/vane_AOA.csv')[start_index])
    q0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bPitchRate.csv')[start_index])
    X0 = np.array([[0],
                   [0],
                   [0],
                   [0]])
    #X0 = np.array([[0.], [0.], [0.], [0.]])
    yout, tout, xout = ml.lsim(sym_state_space(cp), U=U_real, T=t, X0=X0)
    u = []
    alpha = []
    theta = []
    q = []
    for i in range(len(yout)):
        u.append(V0 * (yout[i][0] + 1))
        alpha.append((yout[i][1] + alpha0) * 180 / np.pi)
        theta.append((yout[i][2] + theta0) * 180 / np.pi)
        q.append((yout[i][3]) * 180/ np.pi * V0/ cp.c)
    plt.subplot(2, 3, 1)
    plt.plot(tout, u, label="u, Sim")
    # plt.plot(t, np.multiply(np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index:end_index],
    #                         np.cos(np.genfromtxt('FlightData/vane_AOA.csv')[start_index:end_index])))
    plt.plot(t,np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index:end_index],label = 'data')
    plt.xlabel("Time[s]")
    plt.ylabel('$V_{t}$ [m/s]')
    plt.legend()
    plt.subplot(2, 3, 2)
    plt.plot(tout, alpha, label="alpha, Sim")
    plt.plot(t, np.genfromtxt('FlightData/vane_AOA.csv')[start_index:end_index],label = 'data')
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\alpha$ [deg]')
    plt.legend()
    plt.subplot(2, 3, 3)
    plt.plot(tout, theta, label="theta, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_Pitch.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\theta$ [deg]')
    plt.subplot(2, 3, 4)
    plt.plot(tout, q, label="q, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_bPitchRate.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel('q [deg/s]')
    plt.subplot(2, 3, 5)
    plt.plot(tout, U_real, label="de")
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\delta_{e}$ [deg]')
    plt.show()

#Short Period
def short_period(start_index, end_index, cp):
    t = np.arange(0, 4.1, 0.1)
    #U_real = np.radians(np.genfromtxt('FlightData/delta_e.csv')[35540:35691])
    U_real = np.radians(np.genfromtxt('FlightData/delta_e.csv')[start_index:end_index])
    V0 = np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index]
    theta0 = np.radians(np.genfromtxt('FlightData/Ahrs1_Pitch.csv')[start_index])
    alpha0 = np.radians(np.genfromtxt('FlightData/vane_AOA.csv')[start_index])
    q0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bPitchRate.csv')[start_index])
    X0 = np.array([[0],
                   [0],
                   [0],
                   [0]])
    yout, tout, xout = ml.lsim(sym_state_space(cp), U=U_real, T=t, X0=X0)
    u = []
    alpha = []
    theta = []
    q = []
    for i in range(len(yout)):
        u.append(V0*(yout[i][0]+1))
        alpha.append((yout[i][1]+alpha0) * 180/ np.pi)
        theta.append((yout[i][2]+theta0) * 180/ np.pi)
        q.append((yout[i][3])* 180/ np.pi * V0/ cp.c)
    plt.subplot(2, 3, 1)
    plt.plot(tout, u, label="u, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel('$V_{t}$ [m/s]')
    plt.subplot(2, 3, 2)
    plt.plot(tout, alpha, label="alpha, Sim")
    plt.plot(t, np.genfromtxt('FlightData/vane_AOA.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\alpha$ [deg]')
    plt.subplot(2, 3, 3)
    plt.plot(tout, theta, label="theta, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_Pitch.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\theta$ [deg]')
    plt.subplot(2, 3, 4)
    plt.plot(tout, q, label="q, Sim")
    plt.plot(t,np.genfromtxt('FlightData/Ahrs1_bPitchRate.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel('q [deg/s]')
    plt.subplot(2, 3, 5)
    plt.plot(tout, U_real, label="de")
    plt.legend()
    plt.xlabel("Time[s]")
    plt.ylabel(r'$\delta_{e}$ [deg]')
    plt.show()


# Spiral
def spiral(start_index, end_index, cp):
    t = np.arange(0, 100.1, 0.1)
    #39080:41081

    beta0 = 0.
    V0 = np.radians(np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index])
    phi0 = np.radians(np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index])
    p0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index])
    r0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index])
    X0 = np.array([[0],
                   [0],
                   [0],
                   [0]])
    U = np.transpose(np.array([np.radians(np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index]),
                           np.radians(np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index])]))
    yout, tout, x_out = ml.lsim(asym_state_space(cp), U=U, T =t, X0 = X0)
    beta = []
    phi = []
    p = []
    r = []
    for i in range(len(yout)):
        beta.append((yout[i][0]+beta0) * 180 / np.pi)
        phi.append((yout[i][1]+phi0)* 180 / np.pi)
        p.append((yout[i][2]) * 180 / np.pi * 2 * V0/ cp.b)
        r.append((yout[i][3]) * 180 / np.pi * 2* V0/ cp.b)
    plt.subplot(2, 3, 1)
    plt.plot(tout, beta, label="beta, Sim")
    plt.hlines(beta0,min(t),max(t), label = 'data', colors = 'xkcd:orange')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel(r'$\beta$ [deg]')
    plt.subplot(2, 3, 2)
    plt.plot(tout, phi, label="phi, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index:end_index], label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel(r'$\phi$ [deg]')
    plt.subplot(2, 3, 3)
    #plt.plot(tout, p, label="p, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index:end_index], label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel('p [deg/s]')
    plt.subplot(2, 3, 4)
    plt.plot(tout, r, label="r, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index:end_index], label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel('r [deg/s]')
    plt.subplot(2, 3, 5)
    plt.plot(tout, np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index],label="da")
    plt.plot(tout, np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index], label="dr")
    plt.legend()
    plt.show()
    print(phi)
# Aperiodic Roll
def aperiodic_roll(start_index, end_index, cp):
    #36310:36611
    t = np.arange(0, 60.1, 0.1)
    beta0 = 0.
    V0 = np.radians(np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index])
    phi0 = np.radians(np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index])
    p0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index])
    r0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index])
    X0 = np.array([[0],
                   [0],
                   [0],
                   [0]])
    U = np.transpose(np.array([np.radians(np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index]),
                               np.radians(np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index])]))
    yout, tout, x_out = ml.lsim(asym_state_space(cp), U=U, T=t, X0=X0)
    beta = []
    phi = []
    p = []
    r = []
    for i in range(len(yout)):
        beta.append((yout[i][0]) * 180 / np.pi)
        phi.append((yout[i][1] + phi0) * 180 / np.pi)
        p.append((yout[i][2]))
        r.append((yout[i][3]))
    plt.subplot(2, 3, 1)
    plt.plot(tout, beta, label="beta, Sim")
    plt.hlines(beta0, min(t), max(t), label='data', colors='xkcd:orange')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel(r'$\beta$ [deg]')
    plt.subplot(2, 3, 2)
    plt.plot(tout, phi, label="phi, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index:end_index],label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel(r'$\phi$ [deg]')
    plt.subplot(2, 3, 3)
    plt.plot(tout, p, label="p, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index:end_index], label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel('p [deg/s]')
    plt.subplot(2, 3, 4)
    plt.plot(tout, r, label="r, Sim")
    plt.plot(t, np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index:end_index], label = 'data')
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel('r [deg/s]')
    plt.subplot(2, 3, 5)
    plt.plot(tout, np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index], label="da")
    plt.plot(tout, np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index], label="dr")
    plt.legend()
    plt.show()

# Dutch Roll
def dutch_roll(start_index, end_index, cp):
    #37470:37671
    t = np.arange(0, 30.1, 0.1)
    beta0 = 0.
    V0 = np.radians(np.genfromtxt('FlightData/Dadc1_tas.csv')[start_index])
    phi0 = np.radians(np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index])
    p0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index])
    r0 = np.radians(np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index])
    X0 = np.array([[0],
                   [0],
                   [0],
                   [0]])
    U = np.transpose(np.array([np.radians(np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index]),
                               np.radians(np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index])]))
    yout, tout, x_out = ml.lsim(asym_state_space(cp), U=U, T=t, X0=X0)
    beta = []
    phi = []
    p = []
    r = []
    for i in range(len(yout)):
        beta.append((yout[i][0] + beta0) * 180 / np.pi)
        phi.append((yout[i][1] + phi0) * 180 / np.pi)
        p.append((yout[i][2]) * 180 / np.pi * 2 * V0 / cp.b)
        r.append((yout[i][3]) * 180 / np.pi * 2 * V0 / cp.b)
    plt.subplot(2, 3, 1)
    plt.plot(tout, beta, label="beta, Sim")

    plt.legend()
    plt.subplot(2, 3, 2)
    plt.plot(tout, phi, label="phi, Sim")
    plt.plot(t,np.genfromtxt('FlightData/Ahrs1_Roll.csv')[start_index:end_index])
    plt.legend()
    plt.subplot(2, 3, 3)
    plt.plot(tout, p, label="p, Sim")
    plt.plot(t,np.radians(np.genfromtxt('FlightData/Ahrs1_bRollRate.csv')[start_index:end_index]))
    plt.legend()
    plt.subplot(2, 3, 4)
    plt.plot(tout, r, label="r, Sim")
    plt.plot(t,np.radians(np.genfromtxt('FlightData/Ahrs1_bYawRate.csv')[start_index:end_index]))
    plt.legend()
    plt.subplot(2, 3, 5)
    plt.plot(tout, np.genfromtxt('FlightData/delta_a.csv')[start_index:end_index], label="da")
    plt.plot(tout, np.genfromtxt('FlightData/delta_r.csv')[start_index:end_index], label="dr")
    plt.legend()
    plt.show()
