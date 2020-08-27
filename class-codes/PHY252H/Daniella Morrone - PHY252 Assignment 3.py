#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daniella Morrone
Assignment 3 in PHY252
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


''' Problem 5.92 (b)'''
def f(T):
    h=6.626e-34
    me=9.109e-31
    P=1.013e+5
    k=1.381e-23
    I=13.6*1.602e-19
    A=k/P
    B=(2*np.pi*me*k/(h**2))**(3/2)
    C=-I/k
    return (A*B*(T**(5/2))*np.e**(C/T))-(1/2)

T = np.arange(5000,18000,1)
temp_half=fsolve(f,15000)
line=0*T

plt.plot(T,f(T), color='seagreen')
plt.plot(T,line, color='black', linestyle='--')
plt.plot(temp_half,0, 'ro')
plt.xlabel('T (K)')
plt.ylabel('f(T)')
plt.title('Problem 5.92 (b): f(T) vs. T Plot')


''' Problem 5.92 (d)'''

def x(t):
    h=6.626e-34
    me=9.109e-31
    P=1.013e+5
    I=13.6*1.602e-19
    A=((2*np.pi*me)/(h)**2)
    B=(1/(2*P))*A**(3/2)*I**(5/2)
    C=(4*P/I**(5/2))
    return B*t**(5/2)*np.e**(-1/t)*((1+(C/t**(5/2))*A**(-3/2)*np.e**(1/t))**(1/2)-1)

t = np.arange(0.01,0.3,0.0001)

plt.figure()
plt.plot(t,x(t), color='steelblue')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Problem 5.92 (d): x(t) vs t Plot')

