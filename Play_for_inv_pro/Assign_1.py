import numpy as np
from scipy.optimize import minimize
from seaborn import heatmap
import seaborn as sns
import matplotlib.pyplot as plt
from numba import jit

sns.color_palette("crest", as_cmap=True)


# m_pure = np.zeros((143, 1))
#
# # Formulate real m vector
# for x in range(4, 7):
#     for y in range(1, 9):
#         m_pure[(y * 13 + x)] = (1/6000 - 1/6300) * np.sqrt(2) * 1000

# G_pure = np.zeros((20, 143))
# Figure out the left ray
def G_matrix_form(resolution=1000):  # narrow ray with only one box
    box_dens = int(1000 / resolution)
    num_box_l = 13000 // resolution
    num_box_d = 11000 // resolution
    G_form = np.zeros((20, num_box_l * num_box_d))
    for s in range(10):
        # Figure out the left ray
        G_row_temp = np.zeros((num_box_d, num_box_l))
        seim_place = int(box_dens * (s + 2))
        for i in range(seim_place - 1, -1, -1):
            G_row_temp[(seim_place - 1) - i][i] = 1
        G_row_temp = G_row_temp.flatten()
        G_form[s] = G_row_temp

        # Form the right ray
        G_row_temp = np.zeros((num_box_d, num_box_l))
        for i in range(seim_place, num_box_l, 1):
            G_row_temp[i - seim_place][i] = 1
        G_row_temp = G_row_temp.flatten()
        G_form[s + 10] = G_row_temp
    return G_form


def G_matrix_form_w(resolution=1000):  # narrow ray with wider ray
    box_dens = int(1000 / resolution)
    num_box_l = 13000 // resolution
    num_box_d = 11000 // resolution
    G_form = np.zeros((20, num_box_l * num_box_d))
    for s in range(10):
        # Figure out the left ray
        G_row_temp = np.zeros((num_box_d, num_box_l))
        seim_place = int(box_dens * (s + 2))
        for i in range(seim_place - 1, -1, -1):
            if i - box_dens > 0:
                ray_left_bound = i - box_dens + 1
            else:
                ray_left_bound = 0
            G_row_temp[seim_place - 1 - i][ray_left_bound:i + 1] = 1
        G_row_temp = G_row_temp.flatten()
        G_form[s] = G_row_temp

        # Form the right ray
        G_row_temp = np.zeros((num_box_d, num_box_l))
        for i in range(seim_place, num_box_l, 1):
            if i + box_dens > num_box_l:
                ray_right_bound = num_box_l
            else:
                ray_right_bound = i + box_dens
            G_row_temp[i - seim_place][i:ray_right_bound] = 1
        G_row_temp = G_row_temp.flatten()
        G_form[s + 10] = G_row_temp
    return G_form


# Form pure m
def m_cal(resolution=1000):
    box_dens = int(1000 / resolution)
    num_box_l = 13000 // resolution
    num_box_d = 11000 // resolution
    m_pu = np.zeros((13 * box_dens * 11 * box_dens, 1))
    for x in range(4 * box_dens, 7 * box_dens):
        for y in range(1 * box_dens, 9 * box_dens):
            m_pu[(y * num_box_l + x)] = (1 / 6000 - 1 / 6300) \
                                        * np.sqrt(2) * resolution
    return m_pu


m_pure = m_cal()


# print(G_pure)


def m_est(g_p, m_p):
    t_p = g_p @ m_p
    t_ob = t_p + (np.linalg.norm(t_p) / 20) * np.random.random((20, 1))  # t_pure + noise
    G_tik = lambda epsilon: np.linalg.pinv(g_p.T @ g_p
                            + epsilon ** 2 *
                            np.eye(len(g_p.T))) @ g_p.T
    m_es = lambda epsi: G_tik(epsi) @ t_ob
    diff_t = lambda epsil: (np.linalg.norm((g_p @ m_es(epsil))
                                - t_ob) ** 2) - (20 *
                                ((np.linalg.norm(t_ob) / 20) ** 2))
    eps = minimize(diff_t, 0.01).fun
    print(eps)
    return m_es(eps)


#


G_pure = G_matrix_form(1000)
m_est0 = m_est(G_pure, m_pure)
m_est1 = m_est0.reshape((11, 13))
ax = heatmap(m_est1, linecolor='g')

# m_pure1 = m_pure.reshape((11, 13))
# plt.figure()
# ax2 = heatmap(m_pure1)

t_est = G_pure @ m_est0
plt.figure()
aa = np.linspace(0, 20, 20)
print(aa.shape)
bb = t_est - (G_pure @ m_pure)
print(bb.shape)
cc = np.ones((20, 1)) * (np.linalg.norm(G_pure @ m_pure) / 20)
print(cc.shape)
plt.errorbar(aa, bb, cc)

# ----------higher solution from narrow ray
# G_h = G_matrix_form(200)
# print(G_h.shape)
# m_pure_h = m_cal(200)
# m_est_h = m_est(G_h, m_pure_h)
#
# m_est_h1 = m_est_h.reshape((55, 65))
# plt.figure()
# ax3 = heatmap(m_est_h1)

# -----------higher solution from wide ray
# G_hw = G_matrix_form_w(200)
# m_pure_hw = m_cal(200)
# m_est_hw = m_est(G_hw, m_pure_hw)
#
# m_est_hw1 = m_est_hw.reshape((55, 65))
# plt.figure()
# ax4 = heatmap(m_est_hw1)
# --------------Much higher

# G_hhw = G_matrix_form_w(125)
# m_pure_hhw = m_cal(125)
# m_est_hhw = m_est(G_hhw, m_pure_hhw)
#
# m_est_hhw1 = m_est_hhw.reshape((88, 104))
# plt.figure()
# ax5 = heatmap(m_est_hhw1)

plt.show()
