# -*- coding: utf-8 -*-
"""
Daniella Morrone
PyLab 4 in PHY224

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This is importing the collected data and assigning variables to the dataset
V,I = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/pylab4/PyLab4 Voltage and Current.txt', skiprows=1, unpack=True)

# The following lines define functions g and f as the nonlinear and linear methods
def f(x,a,b):
    return (a * x + b)
def g(x,a,b):
    return (a * x**b)

x = np.linspace(0.1,30,17) 
xlog = np.log(x)
a=0.5882 # Theoretical parameter for physically accurate model fo tungsten
b=3/5 # Theorerical parameter for an ideal black body

# These are the errors of the individual V and I values; used for the nonlinear method
sigmaV = V*0.25/100
sigmaI = I*0.75/100

# These are the V and I values used in the linear method
log_yi = np.log(I)
log_xi = np.log(V)
# These are the errors of the individual V and I values; used for the nonlinear method
log_sigmaV = np.log(sigmaV)
log_sigmaI = np.log(sigmaI)

p_opt, p_cov = curve_fit(f, log_xi , log_yi , (1,0.5882), log_sigmaV , True)
popt, pcov = curve_fit(g, V , I , (1,0.5882), sigmaV , True)


#This is the plot of the linear graph with the data points and the fit line
plt.figure(figsize=(11,6))
plt.scatter(log_xi,log_yi, label="Data Points")
plt.plot(log_xi, f(log_xi,*p_opt),label="Fit Line")
plt.plot(xlog,f(xlog,a,p_opt[1]), label="Theoretical Fit Line for a tungsten")
plt.plot(xlog,f(xlog,b,p_opt[1]), label="Theoretical Fit Line for an ideal blackbody")
plt.errorbar(log_xi,log_yi,yerr=log_sigmaI, fmt=',', ecolor='gray', label ="Error Bars")

plt.legend()
plt.xlim(-0.5,3.7)
plt.title("Linear Method of Current vs. Voltage Graph")
plt.ylabel("Current (mA)")
plt.xlabel("Voltage (V)")

#This is the plot of the nonlinear graph with the data points and the fit line
plt.figure(figsize=(11,6))
plt.scatter(V,I, label="Data Points")
plt.plot(V, g(V,*popt),label="Fit Line")
plt.plot(x,g(x,popt[0],a), label="Theoretical Fit Line for a tungsten")
plt.plot(x,g(x,popt[0],b), label="Theoretical Fit Line for for an ideal blackbody")
plt.errorbar(V,I,yerr=sigmaI, fmt=',', ecolor='gray', label ="Error Bars")
# These error bars are very small

plt.legend()
plt.title("NonLinear Method of Current vs. Voltage Graph")
plt.ylabel("Current (mA)")
plt.xlabel("Voltage (V)")

# This is the plot of the theoretical fit lines on a log x and y-axis graph
plt.figure(figsize=(11,6))
plt.loglog(log_xi, f(log_xi,*p_opt), label="Linear Model Fit Line from the Data")
plt.loglog(V, g(V,*popt), label="Non-Linear Model Fit Line from the Data")

plt.legend()
plt.title("Theoretical Fit Lines on a log y-axis")
plt.ylabel("Log-axis; Current (mA)")
plt.xlabel("Log-axis; Voltage (V)")


''' Modify your program from Section 4 to calculate the standard deviation of the parameters.
What values did you find? Does the value of your fitted exponent fall within the range of the
blackbody values (3/5) with your calculated standard deviation. What about in comparison to
the expected value for tungsten? '''

# a parameter from the linear model
linear_a = p_opt[0]
# b parameter from the nonlinear model
nonlinear_b = popt[1]
# These two parameters represtent the same value in the two different models

# Standard deviation for the a parameter from the linear model
std_l_a = np.sqrt(p_cov[0,0])
# Standard deviation for the b parameter from the nonlinear model
std_nl_b = np.sqrt(pcov[1,1])

'''print("This is linear a (multiply) value:",linear_a,", and its standard deviation:", std_l_a)
print("This is nonlinear b value", nonlinear_b,", and its standard deviation:", std_nl_b)

'''

''' Add a function to your program to calculate χ2red.
What values were computed? How would you interpret your results? '''


def chi2(y_measure,y_predict,y_error):
    return np.sum( (y_measure - y_predict)**2 / y_error**2 )

def chi2red(y_measure, y_predict, y_error, number_of_parameters):
    return chi2(y_measure, y_predict, y_error)/(y_measure.size - number_of_parameters)

# These are the reduced chi squared values for the linear model
# For the ideal blackbody:
red_chi2_lin_bb = chi2red(log_xi,f(xlog,b,p_opt[1]),log_sigmaV,2)
# For the tungsten:
red_chi2_lin_tung = chi2red(log_xi,f(xlog,a,p_opt[1]),log_sigmaV,2)

# These are the reduced chi squared values for the nonlinear model
# For the ideal blackbody:
red_chi2_nonlin_bb = chi2red(V, g(x,popt[0],b),sigmaV,2)
# For the tungsten:
red_chi2_nonlin_tung = chi2red(V, g(x,popt[0],a),sigmaV,2)



