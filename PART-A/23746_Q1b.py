#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 12:05:11 2021

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
    r = np.sqrt(abs(x)**2 + abs(y)**2)
    
    #calculate acceleration in x and y dimension
    ax = -GM*x/r**3
    ay = -GM*y/r**3
    
    #output acceleration 
    return np.array([ax, ay])

def vel(r, v, t):
    "input position and velocity vectors and output velocity vector"
    return v 

def mag(arr):
    "simple function to return magnitude of two dimensional array"
    return np.sqrt(abs(arr[0])**2 + abs(arr[1])**2)


#creating empty lists to assign values to 
x_points = []
vx_points = []
y_points = []
vy_points = []
t_points = []

r = np.array([34.759358,0])  #in AU 
v = np.array([0, 0.185647])  #dimensionless velocity 

def RK_step_r(func, r, v, t, h):
    k1 = h * func(r,v, t)
    k2 = h * func(r + k1/2,v, t+h/2)
    k3 = h * func(r + k2/2,v,  t+h/2)
    k4 = h * func(r + k3, v, t + h)
    
    return (k1 + 2*k2 + 2*k3 + k4)/6

def RK_step_v(func, r, v, t, h):
    k1 = h * func(r,v, t)
    k2 = h * func(r, v + k1/2, t+h/2)
    k3 = h * func(r, v + k2/2,  t+h/2)
    k4 = h * func(r, v + k3, t + h)
    
    return (k1 + 2*k2 + 2*k3 + k4)/6

h_min = 0.0009  #minimum time step 
h_max = 0.02   #maximum time step 
h = 0.01
t = 0 
r_min = 1e-4   #minimum tolerated value for r (used fixed timestep at values close to zero)
t_max = 140
tol = 0.00000000001     #tolerated error 
steps = 0    #count for number of steps 

while t < t_max:    

    r_inter = r + RK_step_r(vel, r, v, t, h) #one step of size h 
    r_single = r_inter + RK_step_r(vel, r, v, t, h)  #step of 2h through two individual steps h 
    r_double = r + RK_step_r(vel, r, v, t, 2*h)  #one step of size 2h 
    
        
    v_inter = v + RK_step_v(acc, r, v, t, h)   #step of h 
    v_single = v_inter + RK_step_v(acc, r, v, t, h)    #step of 2h through two individual steps h 
    v_double = v + RK_step_v(acc, r, v, t, 2*h)    #step of 2h 
    
    #error is difference between step of 2h and two steps of h (error scaled up)
    error = (abs((mag(r_double)) - (mag(r_single)))/abs(mag(r_double)))*10000000000
    #print(error)
    
    if (abs(r_inter[1]) < r_min) | (abs(r_inter[0]) < r_min):  #loop to use a fixed small step size for very small values 
        if h != h_min: 
            h = h_min
    else: 
        
        if (error > tol) and (h > h_min):   
            #decrease/halve step size 
            r = r + RK_step_r(vel, r, v, t, h/2)
            v = v + RK_step_v(acc, r, v, t, h/2)
            h = h/2
            #print('half')
        
        
        
        elif (error < tol) & (h < h_max): 
            #increase/double step size 
            r = r + RK_step_r(vel, r, v, t, 2*h)
            v = v + RK_step_v(acc, r, v, t, 2*h)
            h = h*2
            #print('double')
            
    
        else: 
            #leave step size 
            r = r + RK_step_r(vel, r, v, t, h)
            v = v + RK_step_v(acc, r, v, t, h)
            #print('same')
         
        
    #print(h)
    
    x_points.append(r[0])
    y_points.append(r[1])
    
    vx_points.append(v[0])
    vy_points.append(v[1])
    t_points.append(t)
    
    t = t + h 
    steps = steps + 1

print('Steps per time unit is {:.0f}'.format(steps/t_max))

plt.figure(figsize=(14,10))
plt.plot(x_points, y_points, ls='None', marker='o', ms=0.6, label="Comet's orbit")


plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.title('Orbit of comet around the sun using adaptive time step method', fontsize=16)    
plt.plot(0,0, 'ro', marker='*', ms=10)
plt.legend(markerscale=10., bbox_to_anchor=(0.75, 0.98))
plt.show()
