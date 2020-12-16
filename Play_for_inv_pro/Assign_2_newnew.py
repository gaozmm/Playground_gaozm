import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy import newaxis as na
from scipy.optimize import newton, minimize_scalar
from Jonas_ass2code import *

mars_soil = np.loadtxt('mars_soil.txt').T
mars_soil[1] = -mars_soil[1] + mars_soil[1][0]
mars_peak = find_peaks(mars_soil[1], prominence=550)[0]

plt.plot(mars_soil[0], mars_soil[1], color='r')
plt.title('Mars soil original data')
plt.ylabel('Counts')
plt.xlabel('Velocity (mm/s)')
# plt.scatter(mars_soil[0][mars_peak], mars_soil[1][mars_peak])


peak_num = len(mars_peak)
epsilon = 0.1
alpha = 0.9

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


def dG_gauss_2NP(g, m):  # from jonas
    peaks = int(m.size / 3)
    G = np.empty((np.shape(g)[0], m.size), dtype='float32')

    G[:, :peaks] = (np.exp(-0.5 * (g - m[2 * peaks:, na]) ** 2 / (m[peaks:2 * peaks, na] ** 2)) / (np.sqrt(2 * np.pi) *
                                                                                                   m[peaks:2 * peaks,
                                                                                                   na])).T
    G[:, peaks:2 * peaks] = ((m[:peaks, na] / (np.sqrt(2 * np.pi) * (m[peaks:2 * peaks, na] ** 4))) * np.power(
        (g - m[2 * peaks:, na]), 2)
                             * np.exp(-np.power((g - m[2 * peaks:, na]), 2) / (2 * (m[peaks:2 * peaks, na] ** 2)))
                             - m[:peaks, na] * np.exp(
                -np.power((g - m[2 * peaks:, na]), 2) / (2 * (m[peaks:2 * peaks, na] ** 2))) /
                             (np.sqrt(2 * np.pi) * (m[peaks:2 * peaks, na] ** 2))).T
    G[:, 2 * peaks:] = (m[:peaks, na] * (g - m[2 * peaks:, na]) * np.exp(
        -np.power((g - m[2 * peaks:, na]), 2) / (2 * (m[peaks:2 * peaks, na] ** 2))) / (
                                    np.sqrt(2 * np.pi) * (m[peaks:2 * peaks, na] ** 3))).T
    return G

def F(m_estt, method='gauss'):
    if method == 'gauss':
        de = np.zeros(512)
        for i in range(512):
            for j in range(peak_num):
                de[i] += gauss_cal(mars_soil[0][i], m_estt[j], m_estt[j+2*peak_num], m_estt[j+peak_num])
        return de
    else:
        de = np.zeros(512)
        for i in range(512):
            for j in range(peak_num):
                de[i] += lorentz_cal(mars_soil[0][i], m_estt[j], m_estt[j + 2 * peak_num], m_estt[j + peak_num])
        return de

def object_f(m):
    return np.linalg.norm(mars_soil[1]-F(m))

def grad_f(m):
    G = dG_gauss_2NP(mars_soil[0], m)
    return np.linalg.inv(G.T @ G + (epsilon) ** 2 * np.eye(len(G.T))) @ G.T @ (mars_soil[1] - F(m))

def linesearch(f, df, x0):
    return minimize_scalar(lambda a: f(x0 + a*df), bounds=(-alpha, alpha), method='bounded').x

def m_make():
    m_est = np.zeros(3 * peak_num)
    m_est[:peak_num] = 3000
    m_est[peak_num:2 * peak_num] = 0.25
    m_est[2 * peak_num:3 * peak_num] = mars_soil[0][mars_peak]
    return m_est

m = m_make()

# tol = 1
# max_iter = 100
# for i in range(max_iter):
#     delta_d = mars_soil[1] - F(m)
#     delta_m = grad_f(m)
#     alpha_ = linesearch(object_f, delta_m, m)
#     m += alpha_ * delta_m
#     print(object_f(m))

G_solve = GaussSolver(mars_soil[0], mars_soil[1], m)
result, nor = G_solve.Solve(alpha_bound=0.5)
label_gau = 'Gaussian fit norm' + str(nor[-1]) + 'alpha is 0.5'
plt.plot(mars_soil[0], F(result), label=label_gau)
plt.legend()

plt.figure()
plt.plot(mars_soil[0], mars_soil[1], color='r')
result, nor = G_solve.Solve(alpha_bound=0.3)
label_gau = 'Gaussian fit norm' + str(nor[-1]) + 'alpha is 0.3'
plt.plot(mars_soil[0], F(result), label=label_gau)
plt.legend()

plt.figure()
plt.plot(mars_soil[0], mars_soil[1], color='r')
result, nor = G_solve.Solve(alpha_bound=0.2)
label_gau = 'Gaussian fit norm' + str(nor[-1]) + 'alpha is 0.2'
plt.plot(mars_soil[0], F(result), label=label_gau)
plt.legend()

plt.figure()
plt.plot(mars_soil[0], mars_soil[1], color='r')
result, nor = G_solve.Solve(alpha_bound=0.1)
label_gau = 'Gaussian fit norm' + str(nor[-1]) + 'alpha is 0.1'
plt.plot(mars_soil[0], F(result), label=label_gau)
plt.legend()

plt.figure()
plt.scatter(mars_soil[0][mars_peak], mars_soil[1][mars_peak], s=50, marker='x', c='black')
plt.plot(mars_soil[0], mars_soil[1], color='r')
plt.grid(linestyle='--', c='#D3D3D3')
result, nor = G_solve.Solve(alpha_bound=0.01, epsilon=0.005)
label_gau = 'Gaussian fit norm ' + str(nor)[:6] + ' when alpha is 0.01'
plt.plot(mars_soil[0], F(result), label=label_gau)
plt.legend()

plt.figure()
diff = mars_soil[1] - F(result)
plt.errorbar(mars_soil[0][::8], diff[::8], 0.03*1e4, label='Uncertainty is ' + str(np.std(diff))[:7], ecolor='black',
             c='red', capsize=10)
plt.xlabel('Velocity (mm/s)')
plt.ylabel('Residuals')
plt.legend()


plt.show()