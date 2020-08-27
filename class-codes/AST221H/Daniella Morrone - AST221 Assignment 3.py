#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daniella Morrone
Assignment 3 in AST221

"""

import numpy as np
import matplotlib.pyplot as plt


''' This is part B '''

# This introduces the solar constants for luminosity, radius, and effective temperature
Lo = 3.8395 * (10**26)
Ro = 6.9550826 * (10**8)
To = 57772

# This sets the effective temperature (T) to a certain bound, which is between 2500K and 40000K for the Main Sequence 
T = np.arange(2500,40000,1)
# This creates the equation for the luminosity divided by the sun's luminosity (L_Lo)
# In terms of the effective temperature and solar constants
L_Lo = (((Lo**0.5/(Ro**2))*((T/To)**4))**2)
# Since the luminosity equation is proportional it requires a constant; when taking the log of L, the contsant is added at the end
h = 15
# The constant -10 has been chosen for the plot to fit the ideal y-axis parameters

plt.figure(figsize=(10,8))
# This plots the line for the main sequence
plt.plot(np.log10(T),np.log10(L_Lo)+h, color="black", label = "Main Sequence")
# The following lines change the axes to have the desired limits for the Hertzsprung-Russel Diagram of the Main Sequence
plt.xlim(3.5,4.5)
plt.ylim(-4,6)
plt.gca().invert_xaxis()


''' This is part C '''

# This defines the given equation and inputs the desired values of the radius
sigma = 5.67040040*(10**-8)
Lr = lambda r: 4*np.pi*(r**2)*sigma*(T**4)

# The first radius is 0.1Ro
Lr1 = Lr(0.1*Ro)
plt.plot(np.log10(T),np.log10(Lr1/Lo), color="green", label="Radius of 0.1Ro")
# The second radius is 1Ro
Lr2 = Lr(1*Ro)
plt.plot(np.log10(T),np.log10(Lr2/Lo), color="mediumaquamarine", label="Radius of 1Ro")
# The third radius is 10Ro
Lr3 = Lr(10*Ro)
plt.plot(np.log10(T),np.log10(Lr3/Lo), color="dodgerblue", label="Radius of 10Ro")
# The fourth radius is 100Ro
Lr4 = Lr(100*Ro)
plt.plot(np.log10(T),np.log10(Lr4/Lo), color="mediumblue", label="Radius of 100Ro")


''' This is part D '''

# This adds the point for Betelgeuse onto the plot
plt.scatter(np.log10(3500),np.log10(63000), color="black", label = "Betelgeuse")


''' This is part F '''

# This defines the equation derived from part g) for the lifetime of a star
# Rearranged to solve for luminosity, taking inputs of the lifetime
Lt = lambda t: ((Lo**(3/4))*(10**10)/t)**(4/3)

# The first lifetime is 10Myr divided by the solar luminosity, Lo
Lt1 = np.log10(Lt(10*10**6)/Lo)
plt.hlines(Lt1,3.5,4.5, colors='maroon', label="Lifetime of 10Myr")
# The second lifetime is 100Myr divided by the solar luminosity, Lo
Lt2 = np.log10(Lt(100*10**6)/Lo)
plt.hlines(Lt2,3.5,4.5, colors='red', label="Lifetime of 100Myr")
# The third lifetime is 1000Myr divided by the solar luminosity, Lo
Lt3 = np.log10(Lt(1000*10**6)/Lo)
plt.hlines(Lt3,3.5,4.5, colors='salmon', label="Lifetime of 1000Myr")
# The fourth lifetime is 10000Myr divided by the solar luminosity, Lo
Lt4 = np.log10(Lt(10000*10**6)/Lo)
plt.hlines(Lt4,3.5,4.5, colors='rosybrown', label="Lifetime of 10000Myr")


# The following lines add a legend, axes labels, and a title to the plot
plt.legend(loc='lower left')
plt.xlabel("Temperature (log(K))")
plt.ylabel("Luminosity (log(L/Lo))")
plt.title('Logarithmic Hertzsprung-Russell Diagram of the Main Squence')
