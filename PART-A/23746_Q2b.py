#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 10:03:05 2021

@author: hazelhallett
"""

import numpy as np 
import matplotlib.pyplot as plt

M = 1 #mass of sun in solar mass units 
M_solar = 1.988e30   #mass of sun in kg 
m1 = 1.898e27/M_solar   #mass of Jupiter (body 1)
m2 =  5.683e26/M_solar  #mass of Saturn (body 2)
a1 = 5.204  #semi-major axis in AU Jupiter 
a2 = 9.853 #semi-major axis in AU Saturn 
e1 = 0.049   #Jupiter's eccentricity
e2 = 0.057   #Saturn's eccentricity

def v_aph(a, e):
    P = np.sqrt(a**3)
    v_aph = ((2*np.pi*a)/P)*((1-e)/(1+e))**0.5
    return v_aph
    
#initial conditions 
v1_init = v_aph(a1, e1)
v2_init = v_aph(a2, e2)


def acc(r, v, t):
    "input position as list: r = [x1,x2,y1,y2], input velocity as list: v =[vx1,vx2,vy1,vy2]"
    x1 = r[0]; vx1 = v[0]
    x2 = r[1]; vx2 = v[1]
    y1 = r[2]; vy1 = v[2]
    y2 = r[3]; vy2 = v[3]

    #using mass in solar masses, time in earth years and distances in AU 
    GM = 4*(np.pi)**2
    
    #magnitude of r 
    r1 = np.sqrt(x1**2 + y1**2)
    r2 = np.sqrt(x2**2 + y2**2)
    r21 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    
    
    ax1 = -(GM*x1)/r1**3 + (GM*m2*(x2-x1))/r21**3
    ay1 = -(GM*y1)/r1**3 + (GM*m2*(y2-y1))/r21**3
    
    ax2 = -(GM*x2)/r2**3 - (GM*m1*(x2-x1))/r21**3
    ay2 = -(GM*y2)/r2**3 - (GM*m1*(y2-y1))/r21**3
    
    
    #output acceleration 
    return np.array([ax1, ax2, ay1, ay2])

def vel(r, v, t):
    "input position and velocity vectors and output velocity vector"
    return v 

#function to calculate centre of mass of r coordinate 
def COM_R(r): 
    "r = [x1,x2,y1,y2]"
    x1 = r[0];
    x2 = r[1]; 
    y1 = r[2];
    y2 = r[3]; 

    x_CM = (m1*x1 + m2*x2)
    y_CM = (m1*y1 + m2*y2)
    
    return np.array([x_CM, y_CM, x_CM, y_CM])

#function to calculate centre of mass of v coordinate 
def COM_V(v): 
    "v =[vx1,vx2,vy1,vy2]"
    vx1 = v[0]
    vx2 = v[1]
    vy1 = v[2]
    vy2 = v[3]

    vx_CM = (m1*vx1 + m2*vx2)
    vy_CM = (m1*vy1 + m2*vy2)
    
    return np.array([vx_CM, vy_CM, vx_CM, vy_CM])
    

#creating empty lists to assign values to 

x1_points = []
vx1_points = []
y1_points = []
vy1_points = []

x2_points = []
vx2_points = []
y2_points = []
vy2_points = []

t_points = []


#initial conditions 
r = np.array([a1, a2, 0, 0])  #in AU 
v = np.array([0, 0, v1_init, v2_init])  #dimensionless velocity 

t = 0
h = 0.01
    
while t < 150:
    
    k1 = h * vel(r,v, t)
    k2 = h * vel(r + k1/2,v, t+h/2)
    k3 = h * vel(r + k2/2,v,  t+h/2)
    k4 = h * vel(r + k3, v, t + h)
    
    r = r + (k1 + 2*k2 + 2*k3 + k4)/6
    #position in centre of mass frame 
    r_ = r - (1/(M+m1+m2))*COM_R(r)
    
    x1_points.append(r_[0])
    x2_points.append(r_[1])
    y1_points.append(r_[2])
    y2_points.append(r_[3])

    
    k1 = h * acc(r,v, t)
    k2 = h * acc(r, v + k1/2, t+h/2)
    k3 = h * acc(r, v + k2/2,  t+h/2)
    k4 = h * acc(r, v + k3, t + h)
    
    v = v + (k1 + 2*k2 + 2*k3 + k4)/6
    #velocity in centre of mass frame
    v_ = v - (1/(M+m1+m2))*COM_V(v)

    
    vx1_points.append(v_[0])
    vx2_points.append(v_[1])
    vy1_points.append(v_[2])
    vy2_points.append(v_[3])

    
    t_points.append(t)
    
    t = t + h 
    
    x_COM = r[0]-r_[0]
    y_COM = r[1]-r_[1]
    #print(x_COM)
    #print(y_COM)   

  
#figure formatting
fig, ax = plt.subplots(figsize=(8, 8))
plt.xlim(-10, 10)
plt.ylim(-10,10)
ax.plot(x1_points, y1_points, ls='None', marker='o', ms=0.5, label="Jupiter's orbit")
ax.plot(x2_points, y2_points, ls='None', marker='o', ms=0.5, label="Saturn's orbit")
plt.xlabel('x (AU)', fontsize=20)
plt.ylabel('y (AU)', fontsize=20)
plt.title('Orbit of Jupiter and Saturn around the sun in the centre of mass reference frame', fontsize=16)    


sun, = plt.plot(0,0, 'ro', marker='*', ms=10)
COM, = plt.plot(x_COM, y_COM, marker='*', ms=10)

legend1 = ax.legend(markerscale=10., bbox_to_anchor=(0.75,0.98))
legend2 = ax.legend([COM, sun],['centre of mass ({:.1e},{:.1e})'.format(x_COM, y_COM), 'sun'], loc=3)
ax.add_artist(legend1)
plt.plot(0,0, 'ro', marker='*', ms=10)


plt.show()
