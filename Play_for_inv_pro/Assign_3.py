import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from scipy.interpolate import interp1d
from scipy.stats.contingency import margins
from seaborn import jointplot

gravity_value = np.array([[535, 749, 963, 1177, 1391, 1605, 1819, 2033, 2247, 2461, 2675, 2889],
                          [-15, -24, -31.2, -36.8, -40.8, -42.7, -42.4, -40.9, -37.3, -31.5, -21.8, -12.8]])
gravity_value[1] *= 1e-5
inter_gra_func = interp1d(gravity_value[0], gravity_value[1], fill_value='extrapolate')


G_constant = 6.67430 * 1e-11
r = np.random
r.seed(42)


class MC_MC:
    def __init__(self, num_m=100):
        self.num_m = num_m
        self.c_m = np.linalg.pinv(np.diag((300*np.ones(self.num_m)))** 2)
        self.c_d = np.linalg.pinv(np.diag((1e-5*np.ones(12)))** 2)
        self.d_x = 3420 / num_m
        self.m_0 = np.ones(self.num_m)
        self.inter_dis = np.linspace(0, 3240, self.num_m) + 0.5 * self.d_x
        self.inter_gra = np.clip(inter_gra_func(self.inter_dis), a_min=-1000, a_max=0)

        self.m_0 = (self.inter_gra / (2 * np.pi * G_constant * -1700))
        print(self.m_0)
        grav_ini = self.grav_anomaly(self.m_0, self.num_m)
        plt.plot(self.inter_dis, self.inter_gra)
        plt.plot(gravity_value[0], grav_ini)
        plt.grid(ls='--')
        plt.xlabel('Distance from West edge(m)')
        plt.ylabel('Bouger anomalu (m * s^(-2))')
        plt.figure()

    def get_int_dis(self):
        return self.inter_dis

    def get_m0(self):
        return self.m_0

    def prior_pro(self, m, const=1):
        return const * np.exp(-0.5 * (m-self.m_0).T @ self.c_m @ (m-self.m_0))

    def likelihood_pro(self, m, const=1):
        return const * np.exp(-0.5 * (gravity_value[1] - self.grav_anomaly(m, self.num_m)).T @ self.c_d
                              @ (gravity_value[1] - self.grav_anomaly(m, self.num_m)))

    def walk(self, iteration, step_length):
        #Random walk, with accpetance rate min(1, posterior(perturb)/posterior(current))
        count_acc = 0
        h_bf = self.m_0 + np.random.random(self.num_m)
        h_af = self.m_0 + np.random.random(self.num_m)
        like = np.zeros(iteration)
        post =np.zeros(iteration)
        margin_matrix = np.zeros((self.num_m, iteration)).T
        for p in range(iteration):
            h_bf[0], h_bf[-1] = 0, 0
            h_bf = np.clip(h_bf, a_min=0, a_max=None)
            no_para = np.random.randint(self.num_m)
            like[p] = self.likelihood_pro(h_bf)
            # print('likelihood: ', like[p])
            # print('prior, ', self.prior_pro(h_bf))
            # s_post_pro = -np.log(self.prior_pro(h_bf) * like[p])
            s_post_pro = self.prior_pro(h_bf) * like[p]
            post[p] = np.log10(s_post_pro)
            h_af[no_para] = h_bf[no_para] + step_length * (2*(np.random.random()-0.5))
            # s_post_pro_pert = -np.log(self.prior_pro(h_af) * self.likelihood_pro(h_af))
            s_post_pro_pert = self.prior_pro(h_af) * self.likelihood_pro(h_af)
            if s_post_pro_pert < s_post_pro:
                if s_post_pro_pert/s_post_pro > np.random.random():
                    count_acc += 1
                    h_bf = h_af
                else:
                    h_bf = h_bf
            else:
                count_acc += 1
                h_bf = h_af
            margin_matrix[p] = h_bf

        print('Accept ratio: ', count_acc / iteration)
        print(h_bf)
        return h_bf, like, post, margin_matrix

    def solve(self, number_of_m=False, step_l=0.5, number_iter=5000):
        # if number_of_m is not False:
        #     self.num_m = number_of_m
        #     self.d_x = 3420 / self.num_m
        # num_m = self.num_m
        # m_0 = self.m_0
        # likelihood_pro = self.likelihood_pro
        # prior_pro = self.prior_pro
        # d_x = self.d_x
        # mc = self.walk
        # grav_anomaly = self.grav_anomaly
        # c_m = self.c_m
        # c_d = self.c_d
        m, li, po, mm = self.walk(number_iter, step_l)
        return m, li, po, mm


    @staticmethod
    @njit()
    def grav_anomaly(h_x, num_m, delta=3e-1):
        dx = 3420 / num_m
        delta_g = np.zeros(12)
        for j in range(12):
            x_j = gravity_value[0][j]
            for i in range(num_m):
                x_i = dx * i
                delta_g[j] += G_constant * -1700 * dx * (np.log(((x_i-x_j)**2+h_x[i]**2)/(((x_i-x_j)**2)+delta)))
        return delta_g

num_of_para = 24
MC = MC_MC(num_of_para)

result, likelihood, posterior, mar_matrix = MC.solve(step_l=30, number_iter=60000)
g = MC.grav_anomaly(result, num_of_para) * 1e5
print(g-gravity_value[1]*1e5)
plt.plot(gravity_value[0]*1e-3, g, label='g(m) from MCMC')
plt.grid(ls='--')
plt.scatter(gravity_value[0]*1e-3, gravity_value[1]*1e5, label='Data point')
g_initial = MC.grav_anomaly(MC.get_m0(), num_of_para)*1e5
plt.plot(gravity_value[0]*1e-3, g_initial, label='g(m_0) from initial guess')
plt.xlabel('Distance from West edge(km)')
plt.ylabel('Bouger anomalu (mGal)')
plt.legend()

plt.figure()
plt.scatter(np.linspace(0, len(likelihood), len(likelihood))[::100], -np.log(likelihood)[::100], s=0.3)
plt.xlabel('Iterations')
plt.ylabel('-Log Likelihood')
plt.plot(np.linspace(0, len(likelihood), len(likelihood)), np.ones(len(likelihood))*(0.5*num_of_para), c='grey')
plt.grid(linestyle='--')
plt.figure()
plt.scatter(np.linspace(0, len(posterior), len(posterior))[10000::10], posterior[10000::10], s=0.8, c='r')
plt.ylabel('Posterior in log scale')
plt.xlabel('Iterations discard the burn in phase')
plt.figure()
plt.scatter(np.linspace(0, len(posterior), len(posterior))[10000::10], mar_matrix[10000::10][:, 11],
            s=0.8, c='r')
print(np.corrcoef(np.linspace(0, len(posterior), len(posterior))[10000::10], mar_matrix[10000::10][:, 11]))
plt.xlabel('Iterations discard the burn in phase')
plt.ylabel('Largest value / m')
plt.figure()
fig, axs = plt.subplots(8, 3)
for i in range(8):
    for j in range(3):
        axs[i, j].hist(mar_matrix[10000:][:, i*3+j], histtype='step', bins=30, fill=False, color='black', density=True)
        if i != 7:
            axs[i, j].set_xticklabels(())
plt.figure()
m0, m1 = margins(mar_matrix[10000:])
plt.bar(MC.get_int_dis(), m1.flatten(), edgecolor='black', fill=False)
plt.plot()
plt.xlabel('Distance from West edge(m)')

for i in range(23):
    print('cor', np.corrcoef(mar_matrix[10000:][:, i], mar_matrix[10000:][:, i+1])[0, 1])

# plt.figure()
# fig, axs = plt.subplots(8, 3)
# for i in range(8):
#     for j in range(3):
#         axs[i, j].hist(mar_matrix[10000:][:, i*3+j], histtype='step', bins=30, fill=False, color='black', density=True)
#         if i != 7:
#             axs[i, j].set_xticklabels(())
# for i in range(12):
#     for j in range(2):
#         if i*2+j == 23:
#             break
#         jointplot(x=mar_matrix[10000:][:, i*2+j], y=mar_matrix[10000:][:, i*2+j+1])
#         plt.title('Scatter between'+str(i * 2 + j)+' and '+str(i * 2 + j+1))

# plt.figure()
# result, likelihood = MC.solve(step_l=0.8)
# g = MC.grav_anomaly(result, num_of_para) * 1e5
# plt.plot(np.linspace(0, 3.24, len(g)), g, label='g(m) from MCMC')
# plt.scatter(gravity_value[0]*1e-3, gravity_value[1]*1e5, label='Data point')
# g_initial = MC.grav_anomaly(MC.get_m0(), num_of_para) * 1e5
# plt.plot(np.linspace(0, 3.24, len(g)), g_initial, label='g(m_0) from initial guess')
# plt.xlabel('Distance from West edge(km)')
# plt.ylabel('Bouger anomalu (mGal)')
# plt.legend()
#
# plt.figure()
# result, likelihood = MC.solve(step_l=1.5)
# g = MC.grav_anomaly(result, num_of_para) * 1e5
# plt.plot(np.linspace(0, 3.24, len(g)), g, label='g(m) from MCMC')
# plt.scatter(gravity_value[0]*1e-3, gravity_value[1]*1e5, label='Data point')
# g_initial = MC.grav_anomaly(MC.get_m0(), num_of_para) * 1e5
# plt.plot(np.linspace(0, 3.24, len(g)), g_initial, label='g(m_0) from initial guess')
# plt.xlabel('Distance from West edge(km)')
# plt.ylabel('Bouger anomalu (mGal)')
# plt.legend()

# plt.show()