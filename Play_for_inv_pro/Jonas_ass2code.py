import numpy as np
try:
    import cupy as cp
except:
    print('GPU Computation not available')
import numpy as numpy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from cycler import cycler
from matplotlib.pyplot import rcParams as rc
from matplotlib import font_manager
from numba import njit, cuda
# from cupy import newaxis as na
from scipy.optimize import minimize_scalar
import time

class GaussSolver:
    def __init__(self,x,y,m0,method='CPU'):
        self.method = method
        self.peak_count = int(m0.size / 3)
        self.m = None
        self.norm = None

        if self.method =='CPU':
            try:
                self.x = x.get()
                self.y = y.get()
                self.m0 = m0.get()
                self.G = np.zeros((self.x.size,self.m0.size))

            except:

                self.x = x
                self.y = y
                self.m0 = m0
                self.G = np.zeros((self.x.size, self.m0.size))

        elif self.method =='GPU':
            try:
                self.x = cp.array(x)
                self.y = cp.array(y)
                self.m0 = cp.array(m0)
                self.G = cp.zeros((self.x.size, self.m0.size))

            except:

                self.x = x
                self.y = y
                self.m0 = m0
                self.G = cp.zeros((self.x.size, self.m0.size))

        self.areas = m0[:self.peak_count]
        self.peak_centers = m0[self.peak_count:2*self.peak_count]
        self.widths = m0[self.peak_count*2:]



    @staticmethod
    @njit()
    def dG_CPU(m,G,x):
        peak_count = int(m.size/3)
        areas = np.expand_dims(m[:peak_count],axis=1)
        peak_centers = np.expand_dims(m[peak_count * 2:],axis=1)
        widths = np.expand_dims(m[peak_count:2 * peak_count],axis=1)

        G[:,:peak_count] = (np.exp(-0.5 * (x - peak_centers) ** 2 / (widths ** 2)) / (np.sqrt(2 * np.pi) * widths)).T
        G[:, peak_count:2 * peak_count] = ((areas/(np.sqrt(2*np.pi)*(widths**4)))*np.power((x-peak_centers),2)* np.exp(-np.power((x-peak_centers),2)/(2*(widths**2)))   - areas*np.exp(-np.power((x-peak_centers),2) / (2*(widths**2)))/(np.sqrt(2*np.pi)*(widths**2))).T
        G[:, 2 * peak_count:] = (areas*(x-peak_centers) * np.exp(-np.power((x-peak_centers),2)/(2*(widths**2))) / (np.sqrt(2* np.pi) * (widths**3))).T
        return G

        #
        # G[:,:peak_count] = (widths ** 2 / ((x - peak_centers)**2 + widths**2)).T
        # G[:, peak_count:2 * peak_count] = -(areas * 2 * (x - peak_centers) ** 2 / ((((x-peak_centers)**2/(widths**2))+1)**2 * widths**3)).T
        # G[:, 2 * peak_count:] = (2 * areas * widths **2 * (x-peak_centers) / (((x-peak_centers)**2+widths**2)**2)).T
        # return G

    @staticmethod
    @njit(parallel=True)
    def G_CPU(m,x):
        gauss = lambda a, f, c: a * np.exp(-0.5 * ((x - f) / c) ** 2) / (c * np.sqrt(2 * np.pi))
        # gauss = lambda a, f, c: a * c ** 2 / ((x - f)**2 + c ** 2)  # Lorentz
        result = np.zeros_like(x)
        peak_count = int(m.size/3)
        A = m[:peak_count]
        C = m[peak_count:2 * peak_count]
        F = m[2 * peak_count:]

        for i in range(peak_count):
            result += gauss(A[i], F[i], C[i])
        return result


    def Solve(self,max_iter=1e3,epsilon=0.01,line_search=True,alpha_bound=0.7,tol=1e-4):
        if self.method is 'CPU':
            dG, F = self.dG_CPU, self.G_CPU
        else:
            raise ValueError('Invalid Method')

        if self.method is 'CPU':

            hist = np.linspace(0, 1e4, 10)
            i = 0
            I = np.eye(self.m0.size)
            m = self.m0
            norm = np.zeros(int(max_iter))
            x = self.x
            G = self.G

            def object_f(m):
                return np.linalg.norm(self.y - F(m, x)) ** 2 - 512 * (0.03 * 1e4) ** 2

            def grad_f(m,G):
                return np.linalg.inv(G.T @ G + (epsilon) ** 2 * I) @ G.T @ (self.y - F(m,x))

            def linesearch(f, df, x0):
                return minimize_scalar(lambda a: f(x0 + a * df), bounds=(-alpha_bound, alpha_bound), method='bounded').x


            while np.mean(np.abs(np.diff(hist))) > tol and i < max_iter:

                G = dG(m,G,x)
                delta_d = self.y - F(m,x)
                delta_m = grad_f(m,G)
                if line_search:
                    alpha = linesearch(object_f,delta_m,m)
                else: alpha = alpha_bound
                m += alpha * delta_m
                hist[:-1] = hist[1:]
                hist[-1] = np.linalg.norm(delta_d)
                norm[i] = np.linalg.norm(delta_d)
                i += 1

            return m, norm[norm>0]

