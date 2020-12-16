# Data reading
import numpy as np
import glob
import pandas as pd

ball_d_19 = 0.001 * np.array([19.03, 19.05, 19.10, 19.05])
ball_d_16 = 0.001 * np.array([16.00, 16.00, 16.00, 15.95])
ball_d_11 = 0.001 * np.array([10.95, 11.00, 11.03, 11.00])
height_inc = 0.01 * np.array([25.30, 25.38, 25.29, 25.42])
length_inc = 0.01 * np.array([93.75, 93.69, 93.68, 93.95])
rail_width = 0.001 * np.array([6.19, 6.1, 6.20, 6.15])

# variables without tails were measured by tape, and with tails were by metal ruler
laser_1 = 0.01 * np.array([23.35, 23.38, 23.35, 23.45])
laser_1_2 = 0.01 * np.array([23.37, 23.39, 23.38, 23.36])
laser_2 = 0.01 * np.array([40.15, 40.13, 40.12, 40.15])
laser_2_2 = 0.01 * np.array([40.16, 40.16, 40.14, 40.17])
laser_3 = 0.01 * np.array([55.7, 55.72, 55.73, 55.7])
laser_3_2 = 0.01 * np.array([55.75, 55.78, 55.75, 55.71])
laser_4 = 0.01 * np.array([71.5, 71.49, 71.5, 71.6])
laser_4_2 = 0.01 * np.array([71.52, 71.52, 71.52, 71.52])
laser_5 = 0.01 * np.array([86.25, 86.32, 86.25, 86.31])
laser_5_2 = 0.01 * np.array([86.29, 86.29, 86.28, 86.33])

laser_avg_tape = np.array([np.average(laser_1), np.average(laser_2), np.average(laser_3), np.average(laser_4),
                           np.average(laser_5)])
laser_avg_rule = np.array([np.average(laser_1_2), np.average(laser_2_2), np.average(laser_3_2), np.average(laser_4_2),
                           np.average(laser_5_2)])
laser_avg_all = (laser_avg_rule + laser_avg_tape) / 2

# the angles were measured in degree
tab_ang_e = np.array([14.5, 14.4, 14.9, 14.8])
tab_ang_w = np.array([14.8, 15.9, 15.5, 15.7])
tab_ang_n = np.array([15.7, 15.4, 15.5, 15.2])
tab_ang_s = np.array([14.9, 14.8, 14.9, 14.8])
tab_ang = np.array([0.9, 0.8, 0.5, 0.2])

len_pen = 0.001 * np.array([1890, 1892, 1887.5, 1891.2])

hei_weight_ruler = 0.001 * np.array([18.1, 17.8, 18.0, 18.2])
hei_weight_tri = 0.001 * np.array([16.5, 17.0, 17.5, 17.8])
hei_weight_cali = 0.001 * np.array([18.10, 18.05, 18.05, 18.08])

# Read the data
file_19mm = glob.glob('19mm*')  # Only read the 19mm data
dfs_ball = {str(file): (pd.DataFrame.to_numpy(pd.read_csv(file, header=13))).T for file in file_19mm}
file_14mm = glob.glob('14mm*')
dfs_14 = {str(file): (pd.DataFrame.to_numpy(pd.read_csv(file, header=13))).T for file in file_14mm}
file_9mm = glob.glob('9mm*')
dfs_9 = {str(file): (pd.DataFrame.to_numpy(pd.read_csv(file, header=13))).T for file in file_9mm}
file_pen = glob.glob('*_3c*')
dfs_pen3 = {str(file): (np.loadtxt(file)).T for file in file_pen}
file_pen = glob.glob('*_10c*')
dfs_pen10 = {str(file): (np.loadtxt(file)).T for file in file_pen}
file_pen = glob.glob('*lum.dat')
dfs_pen = {str(file): (np.loadtxt(file)).T for file in file_pen}
file_timing = glob.glob('*Timi*')
dfs_timing = {str(file): (np.loadtxt(file)).T for file in file_timing}

pen_time_3 = np.zeros(4)
count = 0
for key in dfs_pen3:
    pen_time_3[count] = np.average(dfs_pen3[key][1][1:] - dfs_pen3[key][1][:-1])
    count += 1
print(pen_time_3)

pen_time_10 = np.zeros(3)
count = 0
for key in dfs_pen10:
    pen_time_10[count] = np.average(dfs_pen10[key][1][1:] - dfs_pen10[key][1][:-1])
    count += 1

pen_time_r = np.zeros(4)
count = 0
for key in dfs_pen:
    pen_time_r[count] = np.average(dfs_pen[key][1][1:] - dfs_pen[key][1][:-1])
    count += 1

pen_time_t = np.zeros(4)
count = 0
for key in dfs_timing:
    pen_time_t[count] = np.average(dfs_timing[key][1][1:] - dfs_timing[key][1][:-1])
    count += 1