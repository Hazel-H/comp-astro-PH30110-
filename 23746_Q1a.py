#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 11:55:35 2021

"""
import numpy as np
import matplotlib.pyplot as plt


def acc(r, v, t):
    "input position as list: r = [x, y], input velocity as list: v = [vx, vy]"
    x = r[0]
    y = r[1]

    #using mass in solar masses, time in earth years and distances in AU 
    GM = 4*(np.pi)**2
    
    #magnitude of r 
    r = np.sqrt(x**2 + y**2)
    
    ax = -GM*x/r**3  #acceleration in x 
    ay = -GM*y/r**3  #acceleration in y 
    
    #output numpy array of acceleration 
    return np.array([ax, ay])

def vel(r, v, t):
    "input position and velocity vectors and output velocity vector"
    return v 


#creating empty lists to assign values to 
x_points = []
vx_points = []
y_points = []
vy_points = []
t_points = []

#initial conditions
AU = 1.496e11
r_init = 5.2e12/AU  
v_init = 880*(3.154e7/AU) #m/s


r = np.array([r_init,0])  #in AU 
v = np.array([0, v_init])  #dimensionless velocity 

h = 0.01
t = 0 
t_end = 74
steps = 0 

    
while t < t_end:
    
    k1 = h * vel(r,v, t)
    k2 = h * vel(r + k1/2,v, t+h/2)
    k3 = h * vel(r + k2/2,v,  t+h/2)
    k4 = h * vel(r + k3, v, t + h)
    
    
    r = r + (k1 + 2*k2 + 2*k3 + k4)/6
    
    k1 = h * acc(r,v, t)
    k2 = h * acc(r, v + k1/2, t+h/2)
    k3 = h * acc(r, v + k2/2,  t+h/2)
    k4 = h * acc(r, v + k3, t + h)
    
    
    v = v + (k1 + 2*k2 + 2*k3 + k4)/6
    
    x_points.append(r[0])
    y_points.append(r[1])
    
    vx_points.append(v[0])
    vy_points.append(v[1])
    t_points.append(t)
    
    t = t + h 
    steps = steps + 1
    
   
#figure formatting
fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x_points, y_points, ls='None', marker='o', ms=0.5, label="Comet's orbit (h={})".format(h))

plt.xlabel('x (AU)', fontsize=20)
plt.ylabel('y (AU)', fontsize=20)
plt.title('Orbit of comet around the sun', fontsize=20)    
sun, = ax.plot(0,0, 'ro', marker='*', ms=10)
legend1 = ax.legend(markerscale=10., bbox_to_anchor=(0.75,0.98))
legend2 = ax.legend([sun],['sun'], loc=4)
ax.add_artist(legend1)
plt.show()


       
#calculate and print time period of comet 
min_x = (np.min(x_points))
max_x = (np.max(x_points))
a = (max_x - min_x)/2
T = np.sqrt(a**3)
print('Time period of comet is {:.2f}'.format(T))
    
#print number of steps per year
print('Steps per time unit is {:.0f}'.format(steps/t_end))
