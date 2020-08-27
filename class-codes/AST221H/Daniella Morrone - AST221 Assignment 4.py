#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 19:04:10 2019

Daniella Morrone
Assignment 4 in AST221

"""
import numpy as np


''' PROBLEM 1 '''
#1b is shown on paper

h_bar = 1.054571817 *10**-34 #J * s
me = 9.10938356 *10**-31 #kg
mH = 1.6726219 *10**-27 #kg
Mo = 1.989 *10**30 #kg
Ro = 6.9551 *10**8 #m
G = 6.673 *10**-11 #N * m^2 * kg^-2

Z_A = lambda atomic_number, atomic_mass: atomic_number/atomic_mass
Pressure_deg = lambda atomic_number, atomic_mass, density: (((3*np.pi**2)**(2/3))/5) * (h_bar**2 / me) * ( Z_A(atomic_number, atomic_mass) * density/mH)**(5/3)
Density = lambda mass_ratio, radius_ratio: (3 * mass_ratio * Mo)/((4*np.pi)*(radius_ratio*Ro)**3)
Pressure = lambda mass_ratio, radius_ratio: (3 * G * (mass_ratio * Mo)**2)/(8 * np.pi * (radius_ratio * Ro)**4)
Percent_Degeneracy = lambda Degeneracy_Pressure, Pressure: 100*Degeneracy_Pressure/Pressure


''' 1a '''

Z = 1+2 #Atomic number for equal abundance of H (1) and He (2)
A = 1+4 #Atomic mass for equal abundance of H (1) and He (4)
m = 0.1 #Mo
r = 0.1 #Ro

pressure_deg_a = Pressure_deg(Z,A,Density(m,r))
    # 1.6192e+15 N * m^2
pressure_a = Pressure(m,r)
    # 1.3467e+16 N * m^2
percent_deg_a = Percent_Degeneracy(pressure_deg_a,pressure_a)
    # 12.0235 %


''' 1c '''

Z = 2 #Atomic number for He (2), since there is only helium present
A = 4 #Atomic mass for He (4), since there is only helium present
m = 0.1 #Mo
r = 0.01 #Ro

density_prime_c = Density(m,r)
    # 1.4114e+8
pressure_deg_c = Pressure_deg(Z,A,density_prime_c)
    # 1.1949e+20
pressure_prime_c = Pressure(m,r)
    # 1.3467e+20
percent_deg_c = Percent_Degeneracy(pressure_deg_c,pressure_prime_c)
    # 88.7284 %


''' 1d '''

Z = 2 #Atomic number for He (2), since there is only helium present
A = 4 #Atomic mass for He (4), since there is only helium present
m = 0.8 #Mo
r = 0.08 #Ro

density_prime_d = Density(m,r)
    # 2.2052e+6
pressure_deg_d = Pressure_deg(Z,A,density_prime_d)
    # 1.1669e+17
pressure_prime_d = Pressure(m,r)
    # 2.1042e+18
percent_deg_d = Percent_Degeneracy(pressure_deg_d,pressure_prime_d)
    # 5.5455 %


''' PROBLEM 2 '''
# 2a is shown on paper

Msun = 1.9891 *10**30 #kg
Mjupiter = 1.90 *10**27 #kg
Rsun = 6.955 *10**8 #m
Rjupiter = 7.1492 *10**7 #m


''' 2b '''

HourConvert = lambda time_seconds: time_seconds/3600
DayConvert = lambda time_seconds: HourConvert(time_seconds)/24
AUconvert = lambda value: value*6.68458712 *10**-12
Density = lambda mass, radius: (3 * mass)/((4*np.pi)*(radius)**3)
Tidal_Dist = lambda density1, density2, radius1: 2.456 * (density1/density2)**(1/3) * radius1

density_sun = Density(Msun,Rsun)
    # 1411.4861
density_jupiter = Density(Mjupiter,Rjupiter)
    # 1241.3454
tidal_dist = Tidal_Dist(density_sun,density_jupiter,Rsun)
    # 1.7829e+9 m
tidal_dist_AU = AUconvert(tidal_dist) # minimum total distribution in AU
    # 0.01192 AU
tidal_dist_R = tidal_dist/Rsun # minimum total distribution in solar radii
    # 2.5634 Ro
Rsun_AU = AUconvert(Rsun)
    # 0.004649 AU

if Rsun_AU > tidal_dist_AU :
    print('within host star')
else:
    print('not within host star')
# not within the host star


''' 2c '''

Period_Kepler3 = lambda distance, mass: ((4 * np.pi**2 * distance**3)/(G * mass))**0.5

period_sec = Period_Kepler3(tidal_dist,Msun)
    # 41055.4880 s
period_hour = HourConvert(period_sec) # Orbital period in hours
    # 11.4043 hours


''' 2d '''

Radius = 1/AUconvert(1/0.03)
    # 4.4879e+9 m
Period_star_sec = Period_Kepler3(Radius,Msun)
    # 163968.7167 s
Period_star_day = DayConvert(Period_star_sec) # Orbital period of the star in hours
    # 1.8978 days

Orbital_Distance = lambda radius: 2*np.pi*radius

orbit_dist_star = Orbital_Distance(Radius)
    # 2.8199e+10 m
velocity_jupiter = orbit_dist_star/Period_star_sec
    # 1.7198e+5 m/s
velcoity_star = velocity_jupiter * Mjupiter/Msun
    # 164.2716 m/s


'''PROBLEM 3'''

Lo = 3.8395 * 10**33 #erg/s
flux_infrared = 1.41 * 10**4 #erg/s/cm^2
A = 0.343
distanceJ = 7.785472 * 10**11 #m
radiusJ = 6.9911 * 10**7 #m
massJ = 1.90 *10**27 #kg
sigma = 5.6704 * 10**-1 #erg/s/cm


''' 3a '''

Flux = lambda luminosity, distance: luminosity/(4 * np.pi * distance**2)
Luminosity = lambda flux, distance: flux * 4 * np.pi * distance**2
Temperature = lambda flux: (flux/sigma)**(1/4)
Energy_Emitted = lambda flux, distance: sigma * Temperature(flux)**4 * 4 * np.pi * distance**2 

L_infrared = Luminosity(flux_infrared,distanceJ)
E_infrared = Energy_Emitted(flux_infrared,distanceJ)

if L_infrared == E_infrared:
    print('the energy emitted is', E_infrared)
else:
    print('no')
# the energy emitted is 1.0739866202104481e+29
    # 1.0740e+29 erg/s


''' 3b '''

Area_Hemisphere = lambda radius: 2 * np.pi * radius**2
Power = lambda luminosity, distance, radius: Flux(luminosity,distance) * Area_Hemisphere(radius)
Energy_Absorbed = lambda albedo, luminosity, distance, radius: (1-albedo) * Power(luminosity, distance, radius)

energy_absJ = Energy_Absorbed(A,Lo,distanceJ,radiusJ)
    # 1.0170e+25 erg/s


''' 3c '''

Energy_ReEmitted = lambda albedo, luminosity, distance, radius: albedo * Power(luminosity, distance, radius)

energy_totalemit = Power(Lo,distanceJ,radiusJ)
    # 1.5480e+25 erg/s
energy_reemitJ = Energy_ReEmitted(A,Lo,distanceJ,radiusJ)
    # 5.3096e+24 erg/s
excess_energy = np.abs(E_infrared-energy_reemitJ)
    # 1.0739e+29 erg/s


''' 3d '''

massE = 5.972 * 10**24 #kg
radio_energyE = 2.2 * 10**20 #erg/s
mass_earthmantle = 4.01 * 10**24 #kg
mass_jupitercore = 0.14 * massJ #kg

Radioactive_Energy1 = lambda mass_planet1, mass_planet2, radioactive_energy2, mass_core1, mass_core2: (radioactive_energy2 * (mass_planet1/mass_planet2)**4)*(mass_core2/mass_core1)

radioactive_energyJ = Radioactive_Energy1(massJ,massE,radio_energyE, mass_jupitercore, mass_earthmantle)
    # 3.3980e+28 erg/s

''' 3e '''
# Derived on paper

ergConvert = lambda value: value*10**7
Potential = lambda mass, radius: (4 * np.pi * G * mass * Density(mass, radius) * radius**2)/5

grav_potential_RJ = Potential(massJ,radiusJ)
    # 2.0674e+36 J
grav_potential_RJerg = ergConvert(Potential(massJ,radiusJ))
    # 2.0674e+43 erg


''' 3f '''
# equation used explained on paper

mu = 2.2*mH
k = 1.380650424 * 10**-23 #J/K

temperature = mu * mH * (grav_potential_RJ/2)/(3*k)
    # 153609.2329 K
