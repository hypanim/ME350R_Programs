import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

#knowns
l2 = 2.77
w2 = 1
a2 = 0
l3 = 1.95
theta2_start = (140/180)*np.pi
theta2_end = (179/180)*np.pi
points = 1000
theta2_range = np.linspace(theta2_start, theta2_end, points)

#unknowns
theta3_range = []
l1_range = []
w3_range = []
vc_range = []
MA_range = []
a3_range = []
ac_range = []

def positionCalcs(theta2, theta3_guess = 1.1*np.pi, l1_guess = -3):
    def vector_loop(x):
        theta3, l1 = x
        loop1 = l1 - l2*np.cos(theta2) - l3*np.cos(theta3)
        loop2 = l2*np.sin(theta2) + l3*np.sin(theta3)
        return [loop1, loop2]
    
    theta3, l1 = fsolve(vector_loop, [theta3_guess, l1_guess])
    return theta3, l1

#filling in positions
for theta2 in theta2_range:
    temp_theta3, temp_l1 = positionCalcs(theta2)
    theta3_range.append(temp_theta3)
    l1_range.append(temp_l1)

def velocityCalcs(theta2, theta3, w3_guess = -2, vc_guess = -4):
    def vector_loop(x):
        w3, vc = x
        loop1 = vc + l2*w2*np.sin(theta2) + l3*w3*np.sin(theta3)
        loop2 = -l2*w2*np.cos(theta2) - l3*w3*np.cos(theta3)
        return [loop1, loop2]
    
    w3, vc = fsolve(vector_loop, [w3_guess, vc_guess])
    return w3, vc

#filling in velocities
for i in range(len(theta2_range)):
    temp_w3, temp_vc = velocityCalcs(theta2_range[i], theta3_range[i])
    w3_range.append(temp_w3)
    vc_range.append(temp_vc)    

def accelerationCalcs(theta2, theta3, w3, a2, a3_guess = 3, ac_guess = 10):
    def vector_loop(x):
        a3, ac = x
        loop1 = ac + (w2**2)*l2*np.cos(theta2) + a2*l2*np.sin(theta2) + (w3**2)*l3*np.cos(theta3) + a3*l3*np.sin(theta3)
        loop2 = -a2*l2*np.cos(theta2) + (w2**2)*l2*np.sin(theta2) - a3*l3*np.cos(theta3) + (w3**2)*l3*np.sin(theta3)
        return [loop1, loop2]
    
    a3, ac = fsolve(vector_loop, [a3_guess, ac_guess])
    return a3, ac

#filling in acceleration
for i in range(len(theta2_range)):
    temp_a3, temp_ac = accelerationCalcs(theta2_range[i], theta3_range[i], w3_range[i], a2)
    a3_range.append(temp_a3)
    ac_range.append(temp_ac)

#calculating mechanical advantage, magnitude only
for i in range(len(theta2_range)):
    MA_temp = np.abs((w2*l2)/vc_range[i])
    MA_range.append(MA_temp)


#hand calculation values at 148 degrees
theta2_handcalc = 148
vc_handcalc = -4.154
MA_handcalc = 0.667

#plotting
theta2_range_deg = [x*(180/np.pi) for x in theta2_range]
theta3_deg = [x*(180/np.pi) for x in theta3_range]

#plot of sliding velocity of point c vs input angle theta2
plt.axes()
plt.plot(theta2_range_deg, vc_range, label = 'Computational Result')
plt.plot(theta2_handcalc, vc_handcalc, '*', label = 'hand calculation')
plt.xlabel('Theta 2 (degrees)')
plt.ylabel('Sliding Velocity of Point C (m/s)')
plt.title('Sliding Velocity of Point C versus Theta 2')
plt.legend()
plt.show()

#plot of mechanical advantage vs input angle theta2
plt.axes()
plt.plot(theta2_range_deg, MA_range, label = 'Computational Result')
plt.plot(theta2_handcalc, MA_handcalc, '*', label = 'hand calculation')
plt.xlabel('Theta 2 (degrees)')
plt.ylabel('Mechanical Advantage (magnitude)')
plt.title('Mechanical Advantage versus Theta 2')
plt.legend()
plt.show()

#plot of angular velocity of link 3 vs input angle theta2
plt.axes()
plt.plot(theta2_range_deg, w3_range, label = 'Computational Result')
plt.xlabel('Theta 2 (degrees)')
plt.ylabel('Angular velocity of Link 3 (rad/s)')
plt.title('Angular Velocity of Link 3 versus Theta 2')
plt.show()

#plot of sliding acceleration of point c vs input angle theta2
plt.axes()
plt.plot(theta2_range_deg, ac_range, label = 'Computational Result')
plt.xlabel('Theta 2 (degrees)')
plt.ylabel('Acceleration of Point C (m/s^2)')
plt.title('Acceleration of Point C versus Theta 2')
plt.show()