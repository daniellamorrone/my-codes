# -*- coding: utf-8 -*-
"""
Daniella Morrone
PyLab 3 in PHY224

"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This is importing the collected data and assigning variables to the datasets
# This line imports the background data from the decay
samplenumber_bg, countnumber_bg = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/pylab3/RadioactiveDecay_TuesdayOct2_2018_background.txt',skiprows=2,unpack=True)
# This line imports the decay data
sample_number, count_number = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/pylab3/RadioactiveDecay_TuesdayOct2_2018_decay.txt',skiprows=2,unpack=True)

# The following lines define functions g and f as the nonlinear and linear methods
def g(x,a,b):
    return (b*math.e**(a*x))
def f(x,a,c):
    return a*x+c

# Subtract the mean background radiation from the data to get the real count numbers
mean_bg = np.mean(countnumber_bg)
count = count_number - mean_bg

# Calculate the standard error for each data point
num_events = sample_number
num_bg = samplenumber_bg
sigma = np.sqrt(num_events + num_bg)/20

# Convert the count data into rates by dividing divide by the time: 20s
deltat = 20
count_rate = count/deltat
                          

''' Perform the linear regression on (xi; log(yi)) using f as the model functions '''

a = -1/(2.6*60)
yi = count_rate[0]
log_yi = math.log10(count_rate[0])
xi = sample_number*20 - 10
lin = f(xi,a,log_yi)
nonlin = g(xi,a,yi)  

p = np.log(count_rate)
log_sigma = np.absolute(sigma/count_rate)

# This is the plot for the nonlinear data, including the theoretical fit line
plt.figure(figsize=(8,5))

plt.scatter(xi, count_rate, label="Data Points")
p_opt, p_cov = curve_fit(g, xi , count_rate , (a, 0), sigma , True)

plt.plot(xi, g(xi,*p_opt), color="g", label="Theoretical Fit Line")
plt.errorbar(xi,count_rate,yerr=sigma, fmt=',', ecolor='orange', label ="Error Bars")

plt.legend()
plt.title("Data and the Nonlinear Fit of Ba-137m Decay")
plt.xlabel("Time(s)")
plt.ylabel("Number of Geiger Counts")


# This is the plot for the linear data, including the theoretical fit line
plt.figure(figsize=(8,5))

plt.scatter(xi, p, label="Data Points")
p_opt_, p_cov_ = curve_fit(f, xi, p, (a, log_yi), sigma, True)

plt.plot(xi, f(xi,*p_opt_), color="g", label = "Theoretical Fit Line")
plt.errorbar(xi,p,yerr=log_sigma, fmt=',', ecolor='orange', label ="Error Bars")

plt.legend()
plt.title("Data and the Linear Fit of Ba-137m Decay")
plt.xlabel("Time(s)")
plt.ylabel("Number of Geiger Counts")


''' Make a second plot of the same data, but using a logarithmic y axis. One way of doing this is to call 
the pylab.semilogy function. Another way is to call pylab.yscale(`log') after you make the plot. '''

# This is the plot of the theoretical fit lines on a log y-axis graph
plt.figure(figsize=(8,5))
plt.semilogy(xi, g(xi,*p_opt), label="Theoretical Non-Linear Fit Line")
plt.semilogy(xi, f(xi,*p_opt_), label="Theoretical Linear Fit Line")

plt.legend()
plt.title("Theoretical Fit Lines on a log y-axis")
plt.xlabel("Time(s)")
plt.ylabel("Log-axis; Number of Geiger Counts")


''' Which regression method gave a half-life closer to the expected half-life of 2.6 minutes? Can you see the difference on the plots? '''

a_nl,b_nl = p_opt
a_l,c_l = p_opt_

t_half_nl = (math.log(0.5,math.e))/(a_nl*60)
t_half_l = (math.log(0.5,math.e))/(a_l*60)
# The half life determined from the linear function is closer to the real half life value of 2.6 minutes


''' Modify your program from Section 5 to calculate the standard deviation of the half-life.
What values did you find? Does the nominal half-life (2.6 minutes) fall in the range? '''

sigma_hl_nl = np.abs((t_half_nl)*(p_cov[0,0]/a_nl))   # covariance of parameter a divided by parameter a
std_nl_hl = t_half_nl - np.abs(sigma_hl_nl) #standard deviation of nonlinear half life from model

sigma_hl_l = np.abs((t_half_l)*(p_cov_[0,0]/a_l))  #   of parameter a divided by parameter a
std_l_hl = t_half_l - np.abs(sigma_hl_l) #stanard deviation of linear half life from model


''' Add a function to your program to calculate χ2 red. What values were computed? How would you interpret your results? '''

def chi2(y_measure,y_predict,y_error):
    return np.sum( (y_measure - y_predict)**2 / y_error**2 )

def chi2reduced(y_measure, y_predict, y_error, number_of_parameters):
    return chi2(y_measure, y_predict, y_error)/(y_measure.size - number_of_parameters)

# These are the values of reduced chi-squared for the linear and nonlinear methods
RedChi2_Linear = chi2reduced(p, f(xi, *p_opt_), log_sigma, 2)
RedChi2_NonLinear = chi2reduced(count_rate, g(xi, *p_opt), sigma, 2)


