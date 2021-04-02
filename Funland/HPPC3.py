import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

x_strong = np.array([1,2,4,8,16,32,64,128,180, 240,300])
y_time = np.array([21.7199, 11.4898, 5.53075, 4.24043, 4.38414, 3.56064, 2.59841, 1.56359,1.08785,
                   0.837983, 0.773154])
y_speedup = y_time[0] / y_time
x_weak = np.array([1,2,4,8,16,32,64])
y_weak = np.array([0.313521, 0.343963,0.386439,0.550095,1.10603,2.05571,3.19982])
y_w_speed = x_weak / y_weak
def fit_strong(x, p):
    return 1 / ((1 - p) + p / x)

def fit_weak(x, p):
    return 1-p + p*x

popt, pcov = curve_fit(fit_strong, x_strong, y_speedup)
poptw, pcovw = curve_fit(fit_weak, x_weak, y_w_speed)
print(popt)
plt.plot(x_strong, y_speedup, c='b', label='Speed up')
plt.plot(x_strong, fit_strong(x_strong, popt[0]), c='r', ls='--', label='Fit with Amdahl’s law')
plt.grid(linestyle='--')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup t1/tN')
plt.legend()
plt.figure()
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup t1/tN')
plt.grid(linestyle='--')
plt.plot(x_weak, y_w_speed, c='b', label='Speed up')
plt.plot(x_weak, fit_weak(x_weak, poptw[0]), c='r', ls='--', label='Fit with Gustafson’s law')
print(poptw[0])
plt.legend()
plt.show()

