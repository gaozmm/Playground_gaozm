import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy import newaxis as na
from scipy.optimize import newton, minimize_scalar
from numba import njit

mars_soil = np.loadtxt('mars_soil.txt').T
mars_soil[1] = -mars_soil[1] + mars_soil[1][0]
mars_peak = find_peaks(mars_soil[1], prominence=550)[0]

plt.plot(mars_soil[0], mars_soil[1], color='r')
plt.scatter(mars_soil[0][mars_peak], mars_soil[1][mars_peak])
peak_num = len(mars_peak)


def gauss_Dd_Df(z_gf, a_gf, f_gf, c_gf): return a_gf * (z_gf - f_gf) * np.exp(-(z_gf - f_gf) ** 2 / (2 * c_gf ** 2)) \
                                                / (np.sqrt(2 * np.pi) * (c_gf ** 3))


def gauss_Dd_Dc(z_gc, a_gc, f_gc, c_gc): return a_gc * np.exp(-(z_gc - f_gc) ** 2 / (2 * c_gc ** 2)) * (
        -1 + ((z_gc - f_gc) ** 2) / (c_gc ** 2)) / (np.sqrt(2 * np.pi) * c_gc ** 2)


def gauss_Dd_Da(z_ga, a_ga, f_ga, c_ga):
    return np.divide(np.power(np.exp(-(z_ga - f_ga)), 2) * (2 * np.power(c_ga, 2)), (np.sqrt(2 * np.pi) * c_ga))


# def gauss_cal(z_g, a_g, f_g, c_g): return a_g * np.exp((-(z_g - f_g) ** 2) / (2 * c_g ** 2)) / (
#         np.sqrt(2 * np.pi) * c_g)

def gauss_cal(z_g, a_g, f_g, c_g):
    return (a_g / (np.sqrt(2*np.pi) * c_g)) * np.exp(-np.power(z_g - f_g, 2) / (2 * c_g * c_g))


def lorentz_Dd_Da(z_g, a_g, f_g, c_g): return c_g ** 2 / ((z_g - f_g) ** 2 + c_g ** 2)


def lorentz_cal(z_g, a_g, f_g, c_g): return a_g * c_g ** 2 / ((z_g - f_g) ** 2 + c_g ** 2)


def lorentz_Dd_Df(z_g, a_g, f_g, c_g): return a_g * c_g ** 2 * (z_g - f_g) / (((z_g - f_g) ** 2 + c_g ** 2) ** 2)


def lorentz_Dd_Dc(z_g, a_g, f_g, c_g): return (2 * a_g * c_g * ((z_g - f_g) ** 2 + c_g ** 2) - 2 * (
        c_g ** 3) * a_g) / \
                                              (((z_g - f_g) ** 2 + c_g ** 2) ** 2)
from numba import njit

@njit()
def G_par_constrct_gaus(m_vec):
    G_part = np.zeros(3 * peak_num)
    G_part = G_part[np.newaxis, :] + mars_soil[0][:, np.newaxis]
    m_a_new = m_vec[np.newaxis, :peak_num]
    m_f_new = m_vec[np.newaxis, peak_num:2 * peak_num]
    m_c_new = m_vec[np.newaxis, peak_num * 2:]

    G_part[:, :peak_num] = gauss_Dd_Da(G_part[:, :peak_num], m_a_new, m_f_new, m_c_new)
    G_part[:, peak_num:2 * peak_num] = gauss_Dd_Df(G_part[:, peak_num:2 * peak_num], m_a_new, m_f_new, m_c_new)
    G_part[:, peak_num * 2:] = gauss_Dd_Dc(G_part[:, 2 * peak_num:], m_a_new, m_f_new, m_c_new)
    return G_part
    # plt.show()
    # print(G_part)


# ______________________________
def m_make():
    m_est = np.zeros(3 * peak_num)
    m_est[:peak_num] = 3000
    m_est[peak_num:2 * peak_num] = mars_soil[0][mars_peak]
    m_est[2 * peak_num:3 * peak_num] = 0.25
    return m_est
# G_par_constrct_gaus(m_est)

# plt.figure()
# plt.imshow(dG_gauss_2NP(mars_soil[0], m_est))


# ______________________________

def gauss_fit_cal(iteration=1):
    m_est = m_make()
    m_est_new = m_make()
    for iter_time in range(iteration):
        G_part = G_par_constrct_gaus(m_est)
        # m_a_n = m_est[np.newaxis, :peak_num]
        # m_f_n = m_est[np.newaxis, peak_num:2 * peak_num]
        # m_c_n = m_est[np.newaxis, peak_num * 2:]
        def f_cal(m_estt):
            de = np.zeros(512)
            for i in range(512):
                for j in range(peak_num):
                    de[i] += gauss_cal(mars_soil[0][i], m_estt[j], m_estt[j+peak_num], m_estt[j+2*peak_num])

            return de

        d_est = f_cal(m_est)
        # d_est_new = np.zeros(512)
        # for i in range(512):
        #     for j in range(peak_num):
        #         d_est_new[i] += gauss_cal(mars_soil[0][i], m_est_new[j], m_est_new[j + peak_num],
                                          # m_est_new[j + 2 * peak_num])
        # print(d_est)
        delta_d = mars_soil[1] - d_est
        # delta_d_new = mars_soil[1] - d_est_new
        # print(delta_d)
        delta_m = np.linalg.inv(G_part.T @ G_part + (0.001 ** 2) *
                                 np.eye(len(G_part.T))) @ G_part.T @ delta_d
        print('mini is ', np.linalg.norm(delta_d - G_part @ delta_m))

        def diff_cal(alp):
            return np.linalg.norm(mars_soil[1] - f_cal(m_est+alp*delta_m))

        alpha = 0.1
        def line_search(alpha):
            return minimize_scalar(lambda a: diff_cal(a), bounds=(-alpha, alpha), method='bounded').x

        alpha_ = line_search(alpha)
        m_est += alpha_ * delta_m

        # delta_m_new = np.linalg.pinv(G_part_new.T @ G_part_new + (0.1 ** 2) *
        #                          np.eye(len(G_part_new.T))) @ G_part_new.T @ delta_d_new
        # delta_m_new = np.linalg.pinv(G_part_new.T @ G_part_new) @ G_part_new.T @ delta_d_new
        #
        # m_est_new = m_est_new + 0.01*delta_m_new
        #
        # print('No minimize', np.linalg.norm(delta_d_new - G_part_new @ delta_m_new))
        if iter_time % 10 == 1:
            plt.plot(mars_soil[0], d_est)
            plt.plot(mars_soil[0], mars_soil[1])
            plt.show()

    plt.plot(mars_soil[0], d_est)

