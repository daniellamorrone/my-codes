#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daniella Morrone
Assignment 3 in AST222

"""

import numpy as np
import matplotlib.pyplot as plt

# Data for 47 Tuc
RA, Dec, uRA, uDec = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 2/AST222H1S/AST222 Assignment 3/47tuc.txt', unpack=True)
# Data for Keanu
kRA, kDec, kuRA, kuDec = np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 2/AST222H1S/AST222 Assignment 3/keanu.txt', unpack=True)



## CONSTANTS AND FUNCTIONS FOR
        # 2.i
d = 4500 # (pc) distance to 47 Tuc
# Converts speed from arcsec/yr to km/s
Speed_Conv = lambda Mu, Distance: Mu*4.74*Distance

        # 2.iii
# Converts degrees to arcsec
conv_deg_to_as = lambda degree: degree*3600
# Converts arcsec to parsec
conv_as_to_pc = lambda parallax: 1/parallax
# Function for angular distance
ang_dist = lambda Delta_RA, Delta_Dec, Declination: (((Delta_RA * np.cos(np.radians(Declination)))**2 + (Delta_Dec)**2)**(1/2))
# Converts degrees to parsec
conv_deg_to_pc = lambda Angle_Deg, Distance: Distance*np.tan(np.radians(Angle_Deg))

        # 2.iv
G = 4.3009125e-3 # (pc⋅Mo−1⋅(km/s)2) Gravitational constant
# Function for mass from the velocity equation of a cluster
Mass = lambda Velocity_Squared, Radius: Radius*Velocity_Squared/G # Velocity in km/s; Radius in pc
#Gives mass in terms of solar masses

        # 2.v
Mag_47tuc = -9.42 # the absolute magnitude of 47 Tuc
Mag_sun = 4.74 # the sun's bolometric absolute magnitude from textbook
L_sun = 3.8395e+26 # (W) luminosity of the sun from textbook
mass_sun = 1.9891e+30 # (kg) mass of the sun from textbook
# Function for luminosity from the magnitude
Luminosity = lambda Magnitude: 100**((Mag_sun-Magnitude)/5)
# Function to find the mass to light ratio
Gamma = lambda Luminosity, Mass: Mass/Luminosity



''' Problem 2.i '''

vRA = Speed_Conv(uRA*0.001,d) #converts uRA from mas/yr to km/s
vDec = Speed_Conv(uDec*0.001,d) #converts uDec from mas/yr to km/s

# Plot for 47 Tuc positions: Right Ascension vs. Declination in degrees
plt.figure()
plt.title("2.i - 47 Tuc Position Plot")
plt.scatter(RA, Dec, marker='.', c='tomato')
plt.xlabel('Right Ascension (deg)')
plt.ylabel('Declination (deg)')

# Plot for 47 Tuc proper motion: Right Ascension vs. Declination in km/s
plt.figure()
plt.title("2.i - 47 Tuc Proper Motion Plot")
plt.scatter(vRA, vDec, marker='.', c='lightsalmon')
plt.xlabel('Right Ascension (km/s)')
plt.ylabel('Declination (km/s)')


''' Problem 2.ii '''

## Taking the mean of the data gives the mid point thus
        # The centre position in coordinates (Dec, RA)
pos_centre = np.array([np.mean(RA),np.mean(Dec)])
        # = array([  5.98103123, -72.0564269 ])

        # The centre velocity (km/s) in coordinates (vDec, vRA)
v_centre_km = np.array([np.mean(vRA),np.mean(vDec)])
        # = array([117.39522015, -61.53528895])

        # The centre velocity (mas/yr) in coordinates (uDec, uRA)
v_centre_mas = np.array([np.mean(uRA),np.mean(uDec)])
        # = array([ 5.50376091, -2.88491744])


''' Problem 2.iii '''

Dec_centre = Dec-pos_centre[1] # Centres the Dec values at the origin
RA_centre = RA-pos_centre[0] # Centres the RA values at the origin

vDec_centre = vDec-v_centre_km[1] # Centres the vDec values at the origin
vRA_centre = vRA-v_centre_km[0] # Centres the vRA values at the origin

## This is the projected distance of each point
        # It is given from the angular distance and then converted from degrees to parsec 
proj_distance = conv_deg_to_pc(ang_dist(RA_centre,Dec_centre,Dec),d)

## This is the projected velocity of each point
        # It is given from the magnitude of the sum of the two velocity vectors 
proj_velocity = np.sqrt(vDec_centre**2+vRA_centre**2)

# Plot for 47 Tuc velocity and distance: Total Projected Velocity in km/s vs. Projected Distance in pc
plt.figure()
plt.scatter(proj_distance,proj_velocity, marker='.',c='peachpuff')
plt.title("2.iii - 47 Tuc Velocity vs. Distance Plot")
plt.ylabel("Total Projected Velocity (km/s)")
plt.xlabel("Projected Distance (pc)")


''' Problem 2.iv'''

sigma_47tuc_Dec = np.std(vDec_centre)
sigma_47tuc_RA = np.std(vRA_centre)
sigma_47tuc = (sigma_47tuc_RA+sigma_47tuc_Dec)/2

# Using sigma, we get the mean velocity squared
vel_47tuc = 3*sigma_47tuc**2

# Velocity squared is used to determine the mass of each point and we can sum that to get the overall mass in units of solar mass
mass_47tuc = Mass(vel_47tuc,2*np.mean(proj_distance))
        # 574147.9011629332 Mo

''' Problem 2.v '''

# Calculating the luminosity for 47 Tuc from the absolute magnitude, in terms of solar luminosities
L_47tuc = Luminosity(Mag_47tuc)
        # = 461317.5745603791 Lo

# Calculating the mass-to-light ratio for 47 Tuc, in terms of solar mass-to-light ratio
mass_light_47tuc = Gamma(L_47tuc,mass_47tuc)
        # = 1.2445827621245036 γo


''' Problem 2.vi'''

kd = 50000 # (pc) distance to Keanu
Mag_keanu = -12.8 # the absolute magnitude of Keanu


kvRA = Speed_Conv(kuRA*0.001,kd) #converts kuRA from mas/yr to km/s
kvDec = Speed_Conv(kuDec*0.001,kd) #converts kuDec from mas/yr to km/s

## Taking the mean of the data gives the mid point thus
        # The centre position in coordinates (kDec, kRA)
kpos_centre = np.array([np.mean(kRA),np.mean(kDec)])
        # = array([ 78.76169751, -69.02684097])

        # The centre velocity (km/s) in coordinates (kvDec, kvRA)
kv_centre_km = np.array([np.mean(kvRA),np.mean(kvDec)])
        # = array([445.03169942,  58.95057836])

        # The centre velocity (mas/yr) in coordinates (kuDec, kuRA)
kv_centre_mas = np.array([np.mean(kuRA),np.mean(kuDec)])
        # = array([1.87777088, 0.24873662])


kDec_centre = kDec-kpos_centre[1] # Centres the kDec values at the origin
kRA_centre = kRA-kpos_centre[0] # Centres the kRA values at the origin

kvDec_centre = kvDec-kv_centre_km[1] # Centres the kvDec values at the origin
kvRA_centre = kvRA-kv_centre_km[0] # Centres the kvRA values at the origin

## This is the projected distance of each point
        # It is given from the angular distance and then converted from degrees to parsec 
kproj_distance = conv_deg_to_pc(ang_dist(kRA_centre,kDec_centre,kDec),kd)

## This is the projected velocity of each point
        # It is given from the magnitude of the sum of the two velocity vectors 
kproj_velocity = np.sqrt(kvDec_centre**2+kvRA_centre**2)

# Sigma is the average of the standard deviation of both the velocity sets
sigma_keanu_Dec = np.std(kvDec_centre)
sigma_keanu_RA = np.std(kvRA_centre)
sigma_keanu = (sigma_keanu_RA+sigma_keanu_Dec)/2

# Using sigma, we get the mean velocity squared
vel_keanu = 3*sigma_keanu**2

# Velocity squared is used to determine the mass of each point and we can sum that to get the overall mass
mass_keanu = Mass(vel_keanu,2*np.mean(kproj_distance)) 
        # 704835904.9389001 Mo

# Calculating the luminosity for Keanu from the absolute magnitude, in terms of solar luminosities
L_keanu = Luminosity(Mag_keanu) 
        # 10375284.158180127 Lo

# Calculating the mass-to-light ratio for Keanu, in terms of solar mass-to-light ratio
mass_light_keanu = Gamma(L_keanu,mass_keanu)
        # 67.93413020723777 γo

''' Problem 3.i'''

G = 4.3009125e-3 # (pc⋅Mo−1⋅(km/s)2) Gravitational constant
c = 299792 # (km/s) speed of light

Theta_pos = lambda Beta_Coef: (Beta_Coef + (Beta_Coef**2 + 4)**(1/2))/2
Theta_neg = lambda Beta_Coef: (Beta_Coef - (Beta_Coef**2 + 4)**(1/2))/2
Theta_E = lambda Mass, Dd, Ds, Dds: ((Dds/(Ds*Dd))*((4*G*Mass)/(c**2)))**(1/2)
conv_rad_to_mas = lambda Angle_Radians: 206264806.71915*Angle_Radians

def Theta(Beta_Coef, Mass, Dd, Ds, Dds):
    x = np.array([np.abs(Theta_pos(Beta_Coef)), np.abs(Theta_neg(Beta_Coef))])
    return conv_rad_to_mas(min(x)*Theta_E(Mass, Dd, Ds, Dds))

B1 = 1
B2 = 10
B3 = 100
B4 = 1000

M = 1 # Mo
dd = 5e+3 # pc
ds = 10e+3 # pc
dds = 5e+3 # pc

''' Problem 3.ii '''

Th_E = conv_rad_to_mas(Theta_E(M,dd,ds,dds))
        # 1.042040892672401 mas

Th1 = Theta(B1,M,dd,ds,dds)
        # 0.6440166893388252 mas

Th2 = Theta(B2,M,dd,ds,dds)
        # 0.10318238233621209 mas
        
Th3 = Theta(B3,M,dd,ds,dds)
        # 0.01041936709418456 mas

Th4 = Theta(B4,M,dd,ds,dds)
        # 0.0010420398506178874 mas


''' Problem 3.iii '''

Ro = 2.2567e-8 # (pc) radius of sun

Ang_Size = lambda Radius, Distance: 2* np.arctan(Radius/Distance)
delta = conv_rad_to_mas(Ang_Size(Ro,dd)/2)

x = delta/Th_E - Th_E/delta # Beta in terms of Einstein radius
        # -1119.3230749799773

B = x*Th_E # Beta
        # -1166.3804162409524


''' Probelm 3.iv '''


S = 27198952706 # number of stars in mag 20, as found from 'http://www.stargazing.net/david/constel/howmanystars.html'
A = 4*np.pi # surface area in rad^2

A_mas = conv_rad_to_mas(conv_rad_to_mas(A))

star_sep = (A_mas/S)**(1/2)
        # 4433.574139356844


# large, unlikely to see second image

''' Problem 3.v '''

def mu_neg(B,Mass, Dd, Ds, Dds):
    sq = (B**2*(Theta_E(Mass, Dd, Ds, Dds)**2) + 4*(Theta_E(Mass, Dd, Ds, Dds)**2))**(1/2)
    return (B*(Theta_E(Mass, Dd, Ds, Dds))/sq + sq/(B*(Theta_E(Mass, Dd, Ds, Dds))) -2)/4

mu1_neg = mu_neg(B1,M,dd,ds,dds)
mu2_neg = mu_neg(B2,M,dd,ds,dds)
mu3_neg = mu_neg(B3,M,dd,ds,dds)
mu4_neg = mu_neg(B4,M,dd,ds,dds)

