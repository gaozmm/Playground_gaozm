########################################################################################################################
##### Inverse Problems
##### Exercise 1: Vertical Fault -- Part 1
##### Solutions
########################################################################################################################
########################################################################################################################
import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import matrix_rank
from numpy.linalg import inv
from scipy.stats import norm

# (1) Perform discretization, find the matrix G.
# ----------------------------------------------------------------------------------------------------------------------
grav_cons = 6.67*10**(-11)                  # Gravity constant in [m3/ kg*s^2]

# load data, transform x-coord from km to m
# ----------------------------------------------------------------------------------------------------------------------
datafile = open('gravdata.txt')
lines = datafile.readlines()
xk = []
dg = []
for i in lines:
    xk.append(i.split()[0])
    dg.append(i.split('  ')[1])

x = np.asarray(xk, 'float')
x = x*1000
dg = np.asarray(dg, 'float')

# Find matrix G
# ----------------------------------------------------------------------------------------------------------------------
z_base = np.zeros(100)
z_top = np.zeros(100)
G = np.zeros((18, 100))
for j in range(1, 19):
    for i in range(1, 101):
        z_base[i-1] = i*1000
        z_top[i-1] = z_base[i-1]-1000

        G[j-1, i-1] = grav_cons*math.log(((z_base[i-1])**2 + (x[j-1])**2)/((z_top[i-1])**2 + (x[j-1])**2), math.e)


# (4) Find a solution to the above problem, d = Gm, i.e. find _rho.
# Find the inverse of G, G^( g), using Thikonov Regularization:
#  ---------------------------------------------------------------------------------------------------------------------

# Find ||n||^2 = N * sigma^2 with N = Number of measurements = 18.
#  ---------------------------------------------------------------------------------------------------------------------
sigma = 1e-9                                        # Uncertainty of +/- 10^-9 s^-2
n_square = 18*sigma**2                              # noise

# Find epsilon value step 1: coarse logarithmic search to determine order
# of magnitude
#  ---------------------------------------------------------------------------------------------------------------------
epsilonv = np.zeros(40)
diff = np.zeros(40)
m_max = np.zeros(40)
for ep in range(0, 40):
    epsilonv[ep] = 10**(20-ep)                                               # vector containing all epsilon values
    G_g1 = inv(np.dot(G.transpose(), G) + epsilonv[ep]**2 * np.identity(100))
    G_g = np.dot(G_g1, G.transpose())
    m_est = np.dot(G_g, dg)                                                  # calculated densities
    r = np.linalg.norm(dg - np.dot(G, m_est), 2)**2                          # residual
    diff[ep] = r - n_square                                                  # minimizing difference between res and noise -- Should be zero
    m_max[ep] = max(m_est)                                                   # checking maximum density variation

# Plotting epsilon vs.difference
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.subplot(2, 1, 1)
plt.loglog(epsilonv, abs(diff), 'x')
plt.title('logarithmic search for \epsilon')
plt.xlabel('\epsilon')
plt.ylabel('diff')
plt.draw()

# excluding epsilon values resulting in very big density changes
# ----------------------------------------------------------------------------------------------------------------------
kk = np.where(m_max > 1000)
diff = np.delete(diff, kk)
epsilonv = np.delete(epsilonv, kk)

# find epsilon value resulting in smallest difference
# ----------------------------------------------------------------------------------------------------------------------
k = np.argmin(diff)
epsilon = epsilonv[k]

# clearing for second refining search
# ----------------------------------------------------------------------------------------------------------------------
del k, epsilon, epsilonv, diff, m_max, kk

# Find epsilon value step 2: refine search
# ----------------------------------------------------------------------------------------------------------------------
epsilonv = np.zeros(200)
diff = np.zeros(200)
m_max = np.zeros(200)
for ep in range(0, 200):
    epsilonv[ep] = 0.05*(ep+1)*10**(-12)                                         # vector containing all epsilon values
    # G_g1 = inv(np.dot(G.transpose(), G) + epsilonv[ep]**2 * np.identity(100))
    G_g1 = inv(G.transpose() @ G + epsilonv[ep] ** 2 * np.identity(100))

    # G_g = np.dot(G_g1, G.transpose())
    G_g = G_g1 @ G.transpose()

    m_est = G_g @ dg                                                  # calculated densities
    r = np.linalg.norm(dg - np.dot(G, m_est), 2)**2                          # residual
    diff[ep] = r - n_square                                                  # minimizing difference between res and noise --Should be zero
    m_max[ep] = max(m_est)                                                   # checking maximum density variation

# Plotting epsilon vs.difference
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.subplot(2, 1, 2)
plt.plot(epsilonv, abs(diff), 'x')
plt.title('linear search for \epsilon')
plt.xlabel('\epsilon')
plt.ylabel('diff')
plt.draw()


# excluding epsilon values resulting in very big density changes
# ----------------------------------------------------------------------------------------------------------------------
kk = np.where(m_max > 1000)
diff = np.delete(diff, kk)
epsilonv = np.delete(epsilonv, kk)

# find epsilon value resulting in smallest difference
# ----------------------------------------------------------------------------------------------------------------------
k = np.argmin(np.abs(diff))
epsilon = epsilonv[k]

print(epsilon)

# calculate Solution with estimated epsilon
# ----------------------------------------------------------------------------------------------------------------------
epsilon = 6.3*10**(-12)
G_g1 = inv(np.dot(G.transpose(), G) + epsilon**2 * np.identity(100))
G_g = np.dot(G_g1, G.transpose())
m_est = np.dot(G_g, dg)
r = np.linalg.norm(dg - np.dot(G, m_est), 2)**2
diff = r - n_square

# plot solution
# ----------------------------------------------------------------------------------------------------------------------
plt.figure()
plt.plot(m_est, '*')
plt.title('The estimated solution')
plt.xlabel('Depth [km]')
plt.ylabel('m_est [kg/m^3]')
plt.show()




