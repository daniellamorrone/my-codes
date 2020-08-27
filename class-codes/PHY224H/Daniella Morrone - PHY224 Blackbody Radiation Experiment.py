#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 23:20:22 2019

Daniella Morrone
Blackbody Radiation Experiment in PHY224

"""
import numpy as np
import matplotlib.pyplot as plt
import glob2

# Exercise 1: Start

print()
filenames = glob2.glob('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/to_combine/*.txt')  
# This reads the folder with all the files in it and finds all the .txt files within the folder
with open('/Users/daniellamorrone/Desktop/data_file.txt', 'w') as f:
# This opens the output file where all the data in which the individual files will be compiled and calls it f
    for file in filenames:
        with open(file) as infile:
            f.write(infile.read()+'\n')
data = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/data_file.txt')           

print("Exercise 1:")
print("All the data files were combined. Their mean is:", np.mean(data),", and their standard deviation is:",np.std(data))

# Exercise 1: End


# Exercise 2: Start

rm = 2**(1/6)
r = np.linspace(0.1,3,100)
# The linspace numpy function is used here and is given:
# a start value of 0.1 as it cannot include 0 for the equation it is used in
# an end value of 3 as the values past 3 from the eqution it is used in are unnecessary
# the number 100 to say how many numbers will be between 0.1 and 3
V = (rm/r)**12 - 2*(rm/r)**6

fig = plt.figure()
plt.plot(r,V)
plt.xlabel("Radius")
plt.ylabel("Potential (V)")
plt.xlim(0.75,3) # The graph is limited at r=0.75 as the graphed equation has an asymptote at r=1, and limited at r=3 as we do not need values past r=3
plt.ylim(-1.25,1) # The graph is limited at V=-1.25 as the graphed equation does not go below V=-1, and limited at V=1 as the r values continue to get large at an asymptote of r=1, thus are not necessary to include
plt.title("Exercise 2: Lennard-Jones potential")

# Exercise 2: End


# Exercise 3: Start

def MaxBoltz(m,T,v):
    k = 1.38*(10**-23)
    return ((((m/(2*np.pi*k*T))**3)**0.5)*4*np.pi*(v**2)*np.e**(-(m*v**2)/(2*k*T)))
# The function defined as MaxBoltz returns the Maxwell-Boltzmann distribution for particles of imputted mass (m), temperature (T) and velocity (v)

v = np.arange(-40,40,0.1)
# The velocity (v) is set to values bertween -40 and 40 with step size 0.1 between each value
m = 4.652*(10**-23)
T80 = MaxBoltz(m,80,v)
T233 = MaxBoltz(m,233,v)
T298 = MaxBoltz(m,298,v)
# We defined T80, T233, T298 as the MaxBoltz function at varying temperatures but constant mass and velocity 

plt.figure(figsize=(8,5))
plt.plot(v,T80, label="Probability density for Temperature 80K")
plt.plot(v,T233, label="Probability density for Temperature 233K")
plt.plot(v,T298, label="Probability density for Temperature 298K")
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel('Probability Density')
plt.title('Exercise 3: Maxwell-Boltzmann Probability Density for Nitrogen')

# Exercise 3: End


# Exercise 4: Start

print()
data = np.array([[2,5,6,2],
                 [0,1,0,0],
                 [1,1,-1,1],
                 [20,20,1,0],
                 [9,1.6,4.2,2]])

def program(array):
    return array[0]+array[2]
# The function defined as program adds the first and third lines of an inputted array

print("Exercise 4:")
print("This is the sum of lines 1 and 3 from the data:",program(data))

# Exercise 4: End


# Exercise 5: Start

Voltage,Temp,err = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/PYTHON/sample_data.txt',skiprows=3,unpack=True)
# This line loads the data from the given file, skips the first 3 rows and unpacks the columns into 3 variables: Voltage, Temp, err

plt.figure(figsize=(8,5))
plt.scatter(Voltage,Temp)
plt.errorbar(Voltage,Temp,err,fmt='none')
# This line creates the error bars, taking the input values for the Voltage on the x axis, Temp on the y axis and err as the values for the error bars; by setting fmt='none', the line connecting the individual data points is removed
plt.xlabel('Voltage (V)')
plt.ylabel('Temperature (C)')
plt.title('Exercise 5: Simulated data of Temperature vs. Voltage of a black body heated with electric current')

# Exercise 5: End
