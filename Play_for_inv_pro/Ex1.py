import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

gravdata = np.loadtxt('gravdata.txt')
gra_grav = np.reshape(gravdata[:, 1], (18, 1))
x_gra = gravdata[:, 0] * 1000
print(x_gra)
G_cons = 6.67430e-11
delta_z = 1000  # meters
# print(gravdata)
layer_num = 100
delta_rho = np.zeros((100, 1))
g_ele = lambda x, z: 2 * G_cons * z * delta_z / (x ** 2 + z ** 2)
g_ele_dis = lambda x, z: G_cons * np.log(((z + 1000) ** 2 + x ** 2) / (z ** 2 + x ** 2))

sigma = 1e-9  # Uncertainty of +/- 10^-9 s^-2
n_square = 18 * sigma ** 2

G = np.zeros((18, 100))
# for row in range(len(G)):
#     G[row] = np.arange(500, 1000*layer_num + 500, 1000)
#     G[row] = g_ele(x_gra[row], G[row])
# #
for row in range(len(G)):
    G[row] = np.arange(1000, 1000 * (layer_num + 1), 1000)
    G[row] = g_ele_dis(x_gra[row], G[row])


def diff(f, g):
    return np.sum((f - g) ** 2)


# delta_rho = np.linalg.lstsq(G, gra_grav)
epsilon = 1e-10
G_tik = lambda epsilon: np.linalg.inv(G.T @ G + (np.power(epsilon, 2)) * np.eye(100)) @ G.T

delta_rho = lambda epsilon: G_tik(epsilon) @ gra_grav

gra_grav_esi = lambda epsilon: G @ delta_rho(epsilon)

diff_2 = lambda epsilon: abs(np.linalg.norm(gra_grav_esi(epsilon) - gra_grav, 2) ** 2 - n_square)

diff_3 = lambda epsilon: np.sum(gra_grav_esi(epsilon) - n_square)
ep_min = opt.brute(diff_2, ranges=(1e-13, 1e-9))
print(ep_min)
#
#
gra_grav_opt = gra_grav_esi(ep_min.fun)
delta_rho_opt = delta_rho(6.2e-12)
ep = np.linspace(1e-12, 1e-11, 1000)
diffjj = np.linspace(1e-17, 3e-17, 1000)
for i in range(len(ep)):
    diffjj[i] = diff_2(ep[i])
z_la = np.arange(500, 1000 * layer_num + 500, 1000)
fig, ax = plt.subplots(3)
ax[0].plot(x_gra, gra_grav)
ax[0].plot(x_gra, gra_grav_opt)
ax[1].plot(z_la / 1000, delta_rho_opt)

ax[2].plot(ep, diffjj)
plt.show()
