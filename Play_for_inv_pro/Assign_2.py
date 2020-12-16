import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import minimize, curve_fit, minimize_scalar
from Jonas_ass2code import *

mars_soil = np.loadtxt('mars_soil.txt').T
mars_soil[1] = -mars_soil[1] + mars_soil[1][0]
mars_peak = find_peaks(mars_soil[1], prominence=500)[0]
print(mars_soil[0][mars_peak])

plt.plot(mars_soil[0], mars_soil[1], color='r')
plt.scatter(mars_soil[0][mars_peak], mars_soil[1][mars_peak])
# plt.scatter(mars_soil[0], mars_soil[1], linewidths=0.5, marker='x', color='r')
plt.xlabel('Velocity (mm/s)')
plt.ylabel('Counts')

def gauss_Dd_Df(z_gf, a_gf, f_gf, c_gf): return a_gf * (z_gf - f_gf) * np.exp(-(z_gf - f_gf)**2 / (2 * c_gf**2))\
                                                / (np.sqrt(2 * np.pi) * (c_gf ** 3))


def gauss_Dd_Dc(z_gc, a_gc, f_gc, c_gc): return a_gc * np.exp(-(z_gc - f_gc)**2 / (2 * c_gc**2)) * (-1 + ((z_gc - f_gc)**2)/(c_gc**2)) \
                                            / (np.sqrt(2*np.pi) * c_gc**2)


def gauss_Dd_Da(z_ga, a_ga, f_ga, c_ga): return np.exp(-(z_ga - f_ga)**2 / (2 * c_ga**2)) / ((2*np.pi)**0.5 * c_ga)

# print('gaussda is ', gauss_Dd_Da(-12, 2500, -11.9, 0.25))
# exit()

def gauss_cal(z_g, a_g, f_g, c_g): return a_g * np.exp(-(z_g - f_g)**2 / (2 * c_g**2)) / ((2*np.pi)**0.5 * c_g)


def lorentz_Dd_Da(z_g, a_g, f_g, c_g): return c_g ** 2 / ((z_g - f_g)**2 + c_g ** 2)


def lorentz_cal(z_g, a_g, f_g, c_g): return a_g * c_g ** 2 / ((z_g - f_g)**2 + c_g ** 2)


def lorentz_Dd_Df(z_g, a_g, f_g, c_g): return a_g * c_g ** 2 * (z_g - f_g) / (((z_g - f_g)**2 + c_g ** 2)**2)


def lorentz_Dd_Dc(z_g, a_g, f_g, c_g): return (2 * a_g * c_g * ((z_g - f_g)**2 + c_g ** 2) - 2 * (c_g ** 3) * a_g) / \
                                              (((z_g - f_g)**2 + c_g ** 2)**2)


def m_initial_cons():
    m_es = np.zeros(3*len(mars_peak))
    area = 3000
    width = 0.25
    for i in range(len(mars_peak)):
        m_es[3*i] = area
        m_es[3*i+1] = mars_soil[0][mars_peak[i]]
        m_es[3*i+2] = width
    return m_es


m_est = m_initial_cons()


def partial_G_constract(m_vec, method='Gaussian'):
    if method == 'Gaussian':
        Guass_g = True
    else:
        Guass_g = False
    peak_num = len(mars_peak)
    G_part = np.zeros((512, 3 * peak_num))

    if Guass_g is True:
        for no_dat in range(512):
            z_j = mars_soil[0][no_dat]
            for j in range(peak_num):
                G_part[no_dat][3*j] = gauss_Dd_Da(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
                # print(gauss_Dd_Da(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2]))
                G_part[no_dat][3*j+1] = gauss_Dd_Df(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
                G_part[no_dat][3*j+2] = gauss_Dd_Dc(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
    else:
        for no_dat in range(512):
            z_j = mars_soil[0][no_dat]
            for j in range(peak_num):
                G_part[no_dat][3*j] = lorentz_Dd_Da(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
                # print(gauss_Dd_Da(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2]))
                G_part[no_dat][3*j+1] = lorentz_Dd_Df(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
                G_part[no_dat][3*j+2] = lorentz_Dd_Dc(z_j, m_vec[3*j], m_vec[3*j+1], m_vec[3*j+2])
    return G_part

# start iteration

#-------------------------
# def partial_G_cons_sorted(m_vec, method='Gaussian'):
#     peak_num = len(mars_peak)
#     G_part = np.zeros(3 * peak_num)
#     G_part = G_part[np.newaxis, :] + mars_soil[0][:, np.newaxis]
#     print(G_part.shape)
#
#     if method == 'Gaussian':
#         G_part[:][:peak_num] = gauss_Dd_Da()
#
# partial_G_cons_sorted(1)
#-----------------------------
# def diff_cal(G_partial_mat, eps, delta_d, d_obser, d_esti):
#     dm_cal = np.linalg.pinv(G_partial_mat.T @ G_partial_mat + eps ** 2 *
#                           np.eye(len(G_partial_mat.T))) @ G_partial_mat.T @ delta_d
#
#     return d_obser - (d_esti + G_partial_mat@dm_cal)

def object_f(m_est):
    d_d = np.zeros(len(mars_soil[0]))
    for data_point in range(len(mars_soil[0])):
        for peak in range(len(mars_peak)):
            d_d[data_point] += gauss_cal(mars_soil[0][data_point], m_est[3 * peak], m_est[3 * peak + 1],
                                           m_est[3 * peak + 2])
    return np.linalg.norm(mars_soil[1]-d_d)
temp = 1e9
alpha = 0.1
for iter_m in range(400):
    d_est = np.zeros(len(mars_soil[0]))
    for data_point in range(len(mars_soil[0])):
        for peak in range(len(mars_peak)):
            d_est[data_point] += gauss_cal(mars_soil[0][data_point], m_est[3*peak], m_est[3*peak+1], m_est[3*peak+2])

    diff_obs_est = mars_soil[1] - d_est
    print('norm is', np.linalg.norm(diff_obs_est))
    G_par = partial_G_constract(m_vec=m_est)
    if abs(temp - np.linalg.norm(diff_obs_est)) < 1e-1:
        alpha = 0.05
    if abs(temp - np.linalg.norm(diff_obs_est)) < 1e-5:
        alpha = 0.01
    if abs(temp - np.linalg.norm(diff_obs_est)) < 1e-8:
        break

    temp = np.linalg.norm(diff_obs_est)

    def diff_cal(eps):
        dm_cal = np.linalg.pinv(G_par.T @ G_par + (eps ** 2) *
                          np.eye(len(G_par.T))) @ G_par.T @ diff_obs_est
        return np.linalg.norm(mars_soil[1] - (d_est + G_par @ dm_cal)) ** 2 - len(G_par) * 0.03 ** 2
    delta_m = np.linalg.pinv(G_par.T @ G_par + (0.05 ** 2) *
                          np.eye(len(G_par.T))) @ G_par.T @ diff_obs_est
    alp = minimize_scalar(lambda a:object_f(m_est + a*delta_m), bounds=(-alpha, alpha), method='bounded').x
    m_est = m_est + alp * delta_m

plt.imshow(G_par)
plt.figure()
plt.plot(mars_soil[0], mars_soil[1])
plt.plot(mars_soil[0], d_est)

# m_matrix_for_sci = np.zeros((4, 512))+1
# popt, pcov = curve_fit(gauss_cal, mars_soil[0], mars_soil[1])

# plt.plot(mars_soil[0], gauss_cal(mars_soil[0], *popt), 'r-',
#          label='fit: a=%5.3f, f=%5.3f, g=%5.3f' % tuple(popt))
plt.show()