# -*- coding: utf-8 -*-
"""
Daniella Morrone
PyLab 1 in PHY224

"""

import numpy as np
import matplotlib.pyplot as plt

t, pos, err =np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/pylab1/Motion Sensor Position Data Set.txt', skiprows=3, unpack=True)
'''print("This is time:", t)
print("This is position:", pos)
print("This is displacement error:", err)'''

low_displacement=min(pos)
high_displacement=max(pos)
difference = high_displacement-low_displacement
amplitude=difference/2
'''print()
print("This is lower, upper displacement:", low_displacement,",", high_displacement)
print("amplitude is", amplitude, "cm")

print()'''

T=0.70375
mass=0.2004
k=(4*np.pi**2*mass)/T**2

ang_freq=(k/mass)**0.5

'''print("Angular frequency is", ang_freq)  
print()       
print("Q3. k=", k)

print()'''

plt.plot(t, pos)
time = np.array([0,1,2,3,4,5,6])
position = np.array([15,30,1])
plt.title("1.4 Position vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Position (cm)")

t, vel =np.loadtxt('/Users/daniellamorrone/Desktop/school/Second Year/Second Year Semester 1/PHY224H1F/pylab1/Motion Sensor Velocity Data Set.txt', skiprows=3, unpack=True)
'''print("This is time:", t)
print("This is velocity:", vel)'''

plt.figure()
plt.plot(t, vel)

time = np.array([0,1,2,3,4,5,6])
velocity = np.array([15,30,1])
plt.title("1.4 Velocity vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (cm/s)")

plt.figure()

plt.plot(pos, vel)

position = np.array([15,30,1])
velocity = np.array([15,30,1])
plt.title("1.4 Velocity vs. Position Graph")
plt.xlabel("Position (cm)")
plt.ylabel("Velocity (cm/s)")


#This is a trial for loop

# 1.5 programming
position_fake=np.zeros(1001)
velocity_fake=np.zeros(1001)
time_fake=np.linspace(0,10,1001)
dt=0.01
position_fake[0]=amplitude


for i in range(1000):
    position_fake[i+1]=position_fake[i] +dt*velocity_fake[i]
    velocity_fake[i+1]=velocity_fake[i]-dt*ang_freq**2*position_fake[i]

'''print("This is the position for 1.5", position_fake) 
print("This is the velocity for 1.5", velocity_fake) '''

plt.figure()
plt.plot(time_fake,position_fake)     
plt.title("1.5 Position vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Position (cm)")

plt.figure()
plt.plot(time_fake,velocity_fake) 
plt.title("1.5 Velocity vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (cm/s)")

# 1.6) what do we see? we see that the amplitude of the velocity increases over time

plt.figure()
plt.plot(position_fake, velocity_fake)  
plt.title("1.5 Velocity vs. Position Graph")
plt.xlabel("Position (cm)")
plt.ylabel("Velocity (cm/s)")

#1.6) what do we see? A spiral!The displacement in position approaches 0, as does the displacement in velocity
    
#1.6.1 energy

#This is the experimental energy graph
v_ex = vel*0.01 #Converts from cm to m
x_ex = pos*0.01
KE_ex = 0.5*mass*v_ex**2
PE_ex = 0.5*k*x_ex**2
Etot_ex = KE_ex + PE_ex

plt.figure()
plt.title("1.6.1 Kinetic and Potential Energies of the mass-spring system over the recorded time")
plt.plot(t, KE_ex, label = "Kinetic Energy")
plt.plot(t, PE_ex, label = "Potential Energy")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.legend()

plt.figure()
plt.title("1.6.1 Total Energy of the mass-spring system over the recorded time")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.plot(t, Etot_ex)

'''print("1.6.1: This is the maximum experimental energy:",max(Etot_ex),"This is the minimum experiemental energy:",min(Etot_ex))'''

#This is the simulation energy graphs
v = velocity_fake 
x = position_fake
KE = 0.5*mass*v**2
PE = 0.5*k*x**2
Etot = KE + PE

plt.figure()
plt.title("1.6.1 Kinetic and Potential Energies of the simulated mass-spring system over a specific time interval")
plt.plot(time_fake, KE, label = "Kinetic Energy")
plt.plot(time_fake, PE, label = "Potential Energy")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.legend()

plt.figure()
plt.title("1.6.1 Total Energy of the simulated mass-spring system over a specific time interval")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.plot(time_fake, Etot)

'''print("1.6.1: This is the maximum simulated energy:",max(Etot),"This is the minimum simulated energy:",min(Etot))'''

#1.7

fake_position=np.zeros(1001)
fake_velocity=np.zeros(1001)
fake_time=np.linspace(0,10,1001)
dt=0.01
fake_position[0]=amplitude

for i in range(1000):
    fake_position[i+1] = fake_position[i] + dt * fake_velocity[i]
    fake_velocity[i+1] = fake_velocity[i] - dt * (k/mass) * fake_position[i]

'''print("This is the position for 1.7", fake_position) 
print("This is the velocity for 1.7", fake_velocity) '''

plt.figure()
plt.plot(fake_time,fake_position)     
plt.title("1.7 Position vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Position (cm)")

plt.figure()
plt.plot(fake_time,fake_velocity) 
plt.title("1.7 Velocity vs. Time Graph")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (cm/s)")

plt.figure()
plt.plot(fake_position, fake_velocity)  
plt.title("1.7 Velocity vs. Position Graph")
plt.xlabel("Position (cm)")
plt.ylabel("Velocity (cm/s)")


_v = fake_velocity 
_x = fake_position
_KE = 0.5*mass*_v**2
_PE = 0.5*k*_x**2
_Etot = _KE + _PE

plt.figure()
plt.title("1.7 Kinetic and Potential Energies of the simulated mass-spring system over a specific time interval")
plt.plot(fake_time, _KE, label = "Kinetic Energy")
plt.plot(fake_time, _PE, label = "Potential Energy")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.legend()

plt.figure()
plt.title("1.7 Total Energy of the simulated mass-spring system over a specific time interval")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.plot(fake_time, _Etot)

'''print("1.7: This is the maximum simulated energy:",max(_Etot),"This is the minimum simulated energy:",min(_Etot))'''








                  
