import data_read
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import pandas as pd

test = '19mm_w_3.csv'
# To read the data of 19mm rolling balls
# The output is a dict with name: numpy array
# And the array has the form of: first row of time and second row of voltage
tab_ang_e = np.average(data_read.tab_ang_e)
tab_ang_w = np.average(data_read.tab_ang_w)
tab_ang_s = np.average(data_read.tab_ang_s)
tab_ang_n = np.average(data_read.tab_ang_n)


def g_cal_self(acc_b, d_ball, d_rail, inc_angle_deg, h_inc, tab_ang=0):
    r_ball_self = d_ball / 2
    r_rail_self = d_rail / 2
    inc_angle = (inc_angle_deg + tab_ang) * np.pi / 180.0
    return (5 / 7) * (acc_b * h_inc / np.sin(inc_angle)) * (1 / (h_inc - np.sqrt(r_ball_self ** 2 - r_rail_self ** 2) *
                                                                 np.cos(inc_angle)))


def g_cal_troels(acc_b, d_ball, d_rail, inc_angle_deg, tab_ang=0, half=False):
    inc_angle = (inc_angle_deg + tab_ang) * np.pi / 180.0
    if half is True:
        d_rail = d_rail / 2
    return (acc_b / np.sin(inc_angle)) * (1 + 2 * (d_ball ** 2) / (5 * (d_ball ** 2 - d_rail ** 2)))


def peak_find(signal, method='scipy'):
    if method == 'scipy':
        peak_sign = find_peaks(signal, height=1, distance=60)
        return peak_sign
    if method == 'left_peak':
        return 0


def acc_quad_fit(time_peak, s_laser):
    f_quad = lambda t, a, b, c: a * t ** 2 + b * t + c
    val = curve_fit(f_quad, time_peak, s_laser)[0]
    return 2 * val[0]


signal_peak = [x for x in peak_find(data_read.dfs_ball[test][1])[0]]

s_las = data_read.laser_avg_all
acc_fit = acc_quad_fit(data_read.dfs_ball[test][0][signal_peak], s_las)
print(acc_fit)

g_troels = g_cal_troels(acc_fit, np.average(data_read.ball_d_19), np.average(data_read.rail_width), tab_ang_e,
                        half=False)
print(g_troels)
g_troels = g_cal_troels(acc_fit, np.average(data_read.ball_d_19), np.average(data_read.rail_width), tab_ang_w,
                        half=False)
print(g_troels)
g_troels = g_cal_troels(acc_fit, np.average(data_read.ball_d_19), np.average(data_read.rail_width), tab_ang_s,
                        half=False)
print(g_troels)
g_troels = g_cal_troels(acc_fit, np.average(data_read.ball_d_19), np.average(data_read.rail_width), tab_ang_n,
                        half=False)
print(g_troels)


class Pendulum(object):
    def __init__(self, length_pen, length_weight, timing, time_3cm, time_10cm, time_rcm):
        # all the input should be the original arrray
        self.length_pen = length_pen
        self.length_weight = length_weight
        self.time = timing
        self.time3 = time_3cm
        self.time10 = time_10cm
        self.timer = time_rcm
        self.g = 0
        self.uncertainty = 0

    def cal_t(self, option='3'):
        if option == '3':
            time_here = self.time3
        elif option == '10':
            time_here = self.time10
        elif option == 'r':
            time_here = self.timer
        else:
            time_here = self.time
        self.g = np.average(self.length_pen + (self.length_weight / 2)) * (2 * np.pi / np.average(time_here)) ** 2
        self.uncertainty = (2 * np.pi / np.average(time_here)) ** 4 * \
                           ((np.std(self.length_pen) + np.std(self.length_weight) + 0.001) ** 2) + (-2 * np.average(self.length_pen +
                self.length_weight / 2) * (2 * np.pi) ** 2 / (np.average(time_here) ** 3)) ** 2 \
                           * (np.std(time_here) + np.std(self.time)) ** 2

        return self.g, self.uncertainty

pen_me = Pendulum(data_read.len_pen, data_read.hei_weight_cali, data_read.pen_time_t, data_read.pen_time_3,
                  data_read.pen_time_10, data_read.pen_time_r)

print(data_read.pen_time_10)
g_10, un_10 = pen_me.cal_t('10')
print(g_10, un_10)
g_3, un_3 = pen_me.cal_t('3')
print(g_3, un_10)
g_r, un_r = pen_me.cal_t('r')
print(g_r, un_r)


# Calculate by my own calculation
# g_laozi = g_cal_self(acc_fit, np.average(data_read.ball_d_19), np.average(data_read.rail_width),
#                         np.average(data_read.tab_ang_e), np.average(data_read.height_inc))

# print('g is : ', g_troels)
# plt.scatter(dfs['19mm_e_1.csv'][0][signal_peak], dfs['19mm_e_1.csv'][1][signal_peak], marker='x', c='r',linewidths=1)
# plt.plot(dfs['19mm_e_1.csv'][0], dfs['19mm_e_1.csv'][1])
# plt.show()
