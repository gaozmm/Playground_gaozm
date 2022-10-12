import sys
sys.path.append('D:/Labber/Script')
import Labber
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

f = Labber.LogFile('Z:/KianG/20220322_DC_Chip2_Ta/2022/05/Data_0530/NW4_Vbias_gatedown_sweepdown.hdf5')
entry = 6
(Vbias, Vdrop) = f.getTraceXY(entry=entry)
Iac = f.getData('Lockin-2 - Value', entry=entry)  # read data

R_nw = ((abs(Vdrop/1e3)/(abs(Iac)/1e5)))/(7.748092*1e-5)
R_nw = np.clip(R_nw, a_min=None, a_max=3e8)
peaks, _ = find_peaks(R_nw, prominence=30, distance=30)
print(peaks)
peaks = [132, 748]
Vbias *= 1e5
mid_peaks = int((peaks[0] + peaks[1]) / 2)
MARs_arg_right = int((mid_peaks + peaks[1]) / 2)
MARs_arg_left = int((mid_peaks + peaks[0]) / 2)
MARs_arg_right_2 = int((mid_peaks + MARs_arg_right) / 2)
MARs_arg_left_2 = int((mid_peaks + MARs_arg_left) / 2)

plt.grid(linestyle=':')
plt.plot(Vbias, R_nw)
plt.scatter(Vbias[peaks], R_nw[peaks], c='r', s=15)
plt.annotate(text=' ', xy=(Vbias[MARs_arg_right], R_nw[MARs_arg_right]-0.2e8),
             arrowprops=dict(width=5, facecolor='black', shrink=0.05), rotation=180)
plt.annotate(text=' ', xy=(Vbias[MARs_arg_left], R_nw[MARs_arg_left]-0.2e8),
             arrowprops=dict(width=5, facecolor='black', shrink=0.05), rotation=180)
plt.annotate(text=' ', xy=(Vbias[MARs_arg_right_2], R_nw[MARs_arg_right_2]-0.1e8),
             arrowprops=dict(width=5, facecolor='black', shrink=0.05), rotation=180)
plt.annotate(text=' ', xy=(Vbias[MARs_arg_left_2], R_nw[MARs_arg_left_2]-0.1e8),
             arrowprops=dict(width=5, facecolor='black', shrink=0.05), rotation=180)
plt.xlabel('Bias Voltage / $\mu$V')
plt.ylabel('Differential Resistance / $\Omega$')
plt.title('1st and 2nd MARs for NW4_Vbias_gateup_sweepup')
plt.show()