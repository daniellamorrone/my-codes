"""
Daniella Morrone
Assignment 2, Problem 4e in AST222

"""

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0.001, 50, 0.001)
# Defines the distance values between 0.001pc and 50pc with steps of 0.001pc

# These following 3 lines define the functions for the velocity of the black hole, the star cluster, and both combined.
v_bh = lambda x: 131.2*x**(-0.5)
v_sc = lambda x: 131.2*x**(0.05)
v_both = lambda x: v_bh(x) + v_sc(x)

# This creates the actual plot and labels which line is which in a legend
plt.figure(figsize=(8,5))
plt.plot(x,v_bh(x), color='r', label="(i) Circular Velocity of Black Hole")
plt.plot(x,v_sc(x), color='g', label="(ii) Circular Velocity of Star Cluster")
plt.plot(x,v_both(x), color='b', label="(iii) Circular Velocity of Black Hole and Star Cluster")
plt.legend()

# This makes the x axis a logarithmic scale, labels the axes and titles the plot
plt.xscale("log")
plt.xlabel("Distance (logarithmic scale, in pc)")
plt.ylabel("Circular Velocity (km/s)")
plt.title("Circular Velocities vs Distance of a Black Hole, Star Cluster and both combined")


