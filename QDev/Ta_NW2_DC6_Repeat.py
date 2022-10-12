import sys
sys.path.append('D:/Labber/Script')
import Labber
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

f = Labber.LogFile('Z:/KianG/20220711_Ta_DCChip6/2022/07/Data_0720/NW2_currentbias_finer_10reps_sweepup_gatedown.hdf5')
period = 41
R_1 = np.zeros([26, 301])
I_C = np.zeros([10, 26])
for rep in range(10):
    for ent in range(26):
        entry = ent + 41 * rep
        (Ibias, Vdrop) = f.getTraceXY(entry=entry)
        Iac = f.getData('Lockin-2 - Value', entry=entry)  # read data
        R_nw = ((abs(Vdrop/1e3)/(abs(Iac)/1e6)))
        R_nw = np.clip(R_nw, a_min=None, a_max=50e3)
        R_nw[115:185] = 0
        R_nw[0:20], R_nw[280:300] = 0, 0
        peaks, _ = find_peaks(R_nw, prominence=50, distance=100)
        # plt.plot(Ibias, Vdrop)
        # plt.scatter(entry, Ibias[peaks[0]]*1e3)
        # plt.ylabel('$I_C$ / nA')
        I_C[rep][ent] = Ibias[peaks[1]] * 1e3
Ic_std = np.std(I_C, axis=0)
Ic_mean = np.mean(I_C, axis=0)
x = range(26)
plt.rcParams["figure.figsize"] = (10, 5)
plt.plot(x, Ic_mean, 'k-')
plt.fill_between(x, Ic_mean - Ic_std, Ic_mean + Ic_std)
plt.xlabel('$V_{gate}$ / V', fontsize=16)
plt.xlim(0, 25)
plt.xticks([0, 13, 25], [0, -0.78, -1.5], fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('$I_{bias}$ / nA', fontsize=16)
plt.grid(linestyle='--')
plt.title('Average on $I_C$ extraction', fontsize=16)
plt.show()


def heatmap2d(arr: np.ndarray):
    plt.rcParams["figure.figsize"] = (10, 5)
    plt.imshow(arr, cmap='RdBu_r', aspect='auto', interpolation='nearest')
    cb = plt.colorbar(ticks=[0, 25000, 50000])
    cb.set_label(label='Resistance / $k\Omega$', fontsize=14)

#%% The extraction of Ic

# heatmap2d(R_1.T)
# # plt.scatter(range(26), I_C, c='w', s=30)
# plt.plot(range(26), I_C, c='#aabf0a', linewidth=3)
# print(np.arange(0, 240, dtype='float') * -1 / 100)
# plt.xticks([0, 13, 25], [0, -0.78, -1.5])
# plt.yticks([0, 75, 150, 225, 300], [-20, -10, 0, 10, 20])
# plt.xlabel('$V_{gate}$ / V', fontsize=16)
# plt.ylabel('$I_{bias}$ / nA', fontsize=16)
# plt.title('Repeatability $I_C$ extraction', fontsize=16)
#
# plt.show()

#%%