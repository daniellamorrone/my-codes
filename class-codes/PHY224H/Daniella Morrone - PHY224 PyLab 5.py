# -*- coding: utf-8 -*-
"""
Daniella Morrone
PyLab 5 in PHY224

"""
# Import the necessary functions and libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import (poisson,norm)

# Load the data and assign it to variables
sample_number, count_num = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/PyLab 5 files/pylab5/FiestaPlateActivity3sec20min.txt', skiprows=2, unpack=True)
sample_number_bg, count_number_bg = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/PyLab 5 files/pylab5/FiestaBackground3sec20min.txt', skiprows=2, unpack=True)

# The time is represented by the 3 times sample number because the interval that time increases by (delta_time) is 3
time  = sample_number*3
delta_time = 3

# To calculate the actual count number, the mean of the background count data is subtracted from the given count data
mean_count_bg = np.mean(count_number_bg)
count_number = count_num - mean_count_bg

# The rate is determined using the count divided by the time step
# Its uncertainty is calculated as defined in the given equation
rate = count_number / delta_time
uncertainty_rate = np.sqrt(count_number)/delta_time

# This gives the rate with uncertainty
for i in range(len(count_number)):
    a = round(rate[i],1),'+/-',round(uncertainty_rate[i],1)

# This is the uncertainty on the count number
uncertainty_Ns = np.sqrt(count_num + count_number_bg)


''' Plot the histogram of the count (or rate) data.
You should modify the histogram’s appearence, and normalize as needed '''
  
# This plots the histogram of the count data
hist, bin_edges = np.histogram(count_number, bins=40)
plt.hist(count_number, bins=40, density=True, color="lightblue",label="Normalized Plate Activity Data")


''' Add the Poisson probability mass function to qualitatively ﬁt the data.
The most appropriate value for µ can by taking the average value of all of your count data.
This is the Maximum Likelyhood Estimation of the parameter '''

# The norm.fit function gives values for mu and the standard deviation of the Gaussian distribution of the count data
mu_c, std_c = norm.fit(count_number)

# This determines the bin centres and converts them into integers
binCtrs = (bin_edges[1:]+bin_edges[:-1])/2
bCc = np.round(binCtrs)

# This uses the poisson.pmf function to plot the Poisson Probability curve on the histogram
z = poisson.pmf(bCc, mu_c)
plt.plot(bCc,z, color="red",label="Poisson Probability")


''' Add the Gaussian Distribution'''

# This uses the norm.pdf function to plot the Gaussian distribution on the histogram
xmin, xmax = plt.xlim()
x_c = np.linspace(xmin, xmax, 100)
p_c = norm.pdf(x_c, mu_c, std_c)
plt.plot(x_c, p_c, 'k', linewidth=2, color='g',label="Gaussian Distribution")

plt.legend()
plt.xlabel("Count")
plt.ylabel("Probability of measuring the count")
plt.title("Fiesta Plate Count Probability")


''' Are there enough data points for the Poisson and Gaussian functions look the same? '''
# Yes as the plots of the two curves look very similar
# Explained further in the report


''' THE SAME DONE FOR RATE'''
plt.figure()

# This plots the histogram of the rate data
hist, bin_edges = np.histogram(rate, bins=40)
plt.hist(rate, bins=40, density=True, color="lightblue",label="Normalized Plate Activity Data")

# The norm.fit function gives values for mu and the standard deviation of the Gaussian distribution of the rate data
mu_r, std_r = norm.fit(rate)

# This determines the bin centres and converts them into integers
binCtrs = (bin_edges[1:]+bin_edges[:-1])/2
bCr = np.round(binCtrs)

# This uses the poisson.pmf function to plot the Poisson Probability curve on the histogram
z = poisson.pmf(bCr, mu_r)
plt.plot(bCr,z, color="red",label="Poisson Probability")

# This uses the norm.pdf function to plot the Gaussian distribution on the histogram
xmin, xmax = plt.xlim()
x_r = np.linspace(xmin, xmax, 100)
p_r = norm.pdf(x_r, mu_r, std_r)
plt.plot(x_r, p_r, 'k', linewidth=2, color='g',label="Gaussian Distribution")

plt.legend()
plt.xlabel("Rate (counts per second)")
plt.ylabel("Probability of the rate")
plt.title("Fiesta Plate Rate Probability")



''' BONUS '''

# The time is represented by the 3 times background sample number because the interval that time increases by (delta_time_bg) is 3
time_bg = sample_number_bg*3
delta_time_bg = 3

# The background rate is determined using the background count divided by the time step
# Its uncertainty is calculated as defined in the given equation
rate_bg = count_number_bg / delta_time_bg
uncertainty_rate_bg = np.sqrt(count_number_bg)/delta_time_bg

# This gives the background rate with uncertainty
for i in range(len(count_number_bg)):
    b = round(rate_bg[i],1),'+/-',round(uncertainty_rate_bg[i],1)

''' Plot the histogram of the count (or rate) data.
You should modify the histogram’s appearence, and normalize as needed '''     

plt.figure()

# This plots the histogram of the background count data
hist_bg, bin_edges_bg = np.histogram(count_number_bg, bins=10)
plt.hist(count_number_bg, bins=10, density=True, color="lightblue", label="Normalized Plate Background-Activity Data")


''' Add the Poisson probability mass function to qualitatively ﬁt the data.
The most appropriate value for µ can by taking the average value of all of your count data.
This is the Maximum Likelyhood Estimation of the parameter ''' 

# The norm.fit function gives values for mu and the standard deviation of the Gaussian distribution of the background count data
mu_cbg, std_cbg = norm.fit(count_number_bg)

# This determines the bin centres and converts them into integers
binCtrs = (bin_edges_bg[1:]+bin_edges_bg[:-1])/2
bCc_bg = np.round(binCtrs)

# This uses the poisson.pmf function to plot the Poisson Probability curve on the histogram
z_bg = poisson.pmf(bCc_bg, mu_cbg)
plt.plot(bCc_bg, z_bg, color="red",label="Poisson Probability")


''' Add the Guassian Distribution '''

# This uses the norm.pdf function to plot the Gaussian distribution on the histogram
xmin, xmax = plt.xlim()
x_cbg = np.linspace(xmin, xmax, 100)
p_cbg = norm.pdf(x_cbg, mu_cbg, std_cbg)
plt.plot(x_cbg, p_cbg, 'k', linewidth=2, color="g",label="Gaussian Distribution")

plt.legend()
plt.xlabel("Background Count")
plt.ylabel("Probability of measuring the count")
plt.title("Fiesta Plate Background Count Probability")


''' Are there enough data points for the Poisson and Gaussian functions look the same? '''
# No as the plot of the poisson tends towards the bars of the histogram and the gaussian is more of a curve following the trend of the data
# Explained further in the report


''' THE SAME DONE FOR BACKGROUND RATE'''
plt.figure()

# This plots the histogram of the background rate data
hist, bin_edges = np.histogram(rate_bg, bins=10)
plt.hist(rate_bg, bins=10, density=True, color="lightblue",label="Normalized Plate Activity Data")

# The norm.fit function gives values for mu and the standard deviation of the Gaussian distribution of the background rate data
mu_rbg, std_rbg = norm.fit(rate_bg)

# This determines the bin centres and converts them into integers
binCtrs = (bin_edges[1:]+bin_edges[:-1])/2
bCr_bg = np.round(binCtrs)

# This uses the poisson.pmf function to plot the Poisson Probability curve on the histogram
z = poisson.pmf(bCr_bg, mu_rbg)
plt.plot(bCr_bg,z, color="red",label="Poisson Probability")

# This uses the norm.pdf function to plot the Gaussian distribution on the histogram
xmin, xmax = plt.xlim()
x_rbg = np.linspace(xmin, xmax, 100)
p_rbg = norm.pdf(x_rbg, mu_rbg, std_rbg)
plt.plot(x_rbg, p_rbg, 'k', linewidth=2, color='g',label="Gaussian Distribution")

plt.legend()
plt.xlabel("Background Rate (counts per second)")
plt.ylabel("Probability of the rate")
plt.title("Fiesta Plate Background Rate Probability")




