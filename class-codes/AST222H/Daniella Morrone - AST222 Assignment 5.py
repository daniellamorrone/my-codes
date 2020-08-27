"""
Daniella Morrone
Assignment 5 in AST222

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

PropDist,VelAway = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 2/AST222H1S/MorroneDaniella.txt', unpack=True)

plt.figure(figsize=(8,5))
plt.scatter(PropDist,VelAway, color='steelblue', marker='+', label='Galaxy Data Points')

vel = lambda Ho, d: Ho*d

Ho,p_cov = curve_fit(vel, PropDist, VelAway)
plt.plot(PropDist,vel(*Ho,PropDist), color='lightcoral', linewidth=0.5, label='Linear Fit')

plt.title('Velocity vs. Proper Distance of Galaxies in the Daniella Morrone Universe')
plt.xlabel('Proper Distance (Mpc)')
plt.ylabel('Velocity moving Away (km/s)')
plt.legend()

Conv_Mpc_km = lambda Value_in_Mpc: Value_in_Mpc*3.086e+19
Conv_s_Gyr = lambda Value_in_s: (((((Value_in_s/60)/60)/24)/365)/1000000000)

inverse_Ho_in_s = Conv_Mpc_km(1/Ho)
Age = Conv_s_Gyr(inverse_Ho_in_s)

G = 6.6742867e-11
CriticalDensity = (3*(1/inverse_Ho_in_s)**2)/(8*np.pi*G)
Dense = 2.56e-27

plt.text(570, 100000, r'The Hubble parameter is', fontsize=10)
plt.text(630, 82500, np.round(*Ho,1),fontsize=10)
plt.text(705, 82500, r'km/s/Mpc.', fontsize=10)

plt.text(515, 58000, r'The age of the Universe is', fontsize=10)
plt.text(843, 58000, np.round(*Age,1), fontsize=10)
plt.text(885, 58000, r'Gyr.', fontsize=10)

plt.text(397, 33000, r'The critical denisty of the Universe is', fontsize=10)
plt.text(860, 33000, '{:.2g}'.format(*CriticalDensity), fontsize=10)
plt.text(955, 33000, r'kg/m^3.', fontsize=10)
plt.text(525, 15500, r'Therefore, the universe is', fontsize=10)
plt.text(848, 15500,
         r'open.' if np.round(CriticalDensity,26) > np.round(Dense,26) else r'closed.' if np.round(CriticalDensity,26) < np.round(Dense,26) else r'flat.',
         fontsize=10)



