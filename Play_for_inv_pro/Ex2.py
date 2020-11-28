import numpy as np
import matplotlib.pyplot as plt

gravdata = np.loadtxt('gravdata.txt')
gra_grav = np.reshape(gravdata[:,1], (18,1))
x_gra = gravdata[:,0] * 1000
print(x_gra)
G_cons = 6.67430e-11



g_ele_dis = lambda x, z:G_cons * np.log(((z+1000)**2 + x**2)/(z**2 + x**2))
