'''
Daniella Morrone
Assignment 1 in PHY252
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


R = 8.3145
P = 101325
T = np.array([100,200,300,400,500,600])
B = np.array([-160,-35,-4.2,9,16.9,21.3])
sigma = np.array([0,0,0,0,0,0])

plt.scatter(T,B, color="r",label="Data Points from (a)")


B_T = lambda Temperature, a, b : b - a/(R*Temperature)
Temp = np.arange(100,600)
p_opt, p_cov = curve_fit(B_T, T , B, (1.86e+5,64), sigma, True)

plt.plot(Temp,B_T(Temp,*p_opt), color="b",label="Curve Fit, Trial and Error")
plt.legend()
plt.ylabel("B(T) (cm^3/mol)")
plt.xlabel("Temperature (K)")

a=207862.5
b=90

plt.plot(Temp,B_T(Temp,a,b), color="g",label="Curve Fit, Calculated")
plt.legend()
plt.ylabel("B(T) (cm^3/mol)")
plt.xlabel("Temperature (K)")
