import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

def positionCalcs(l2, l3, l4, theta2, theta3_guess = 1, l1_guess = 1):
    def vector_loop(x):
        theta3, l1 = x
        loop1 = l2*np.cos(theta2) + l3*np.cos(theta3) - l1
        loop2 = l2*np.sin(theta2) + l3*np.sin(theta3) -l4
        return[loop1, loop2]

    theta3, l1 = fsolve(vector_loop, [theta3_guess, l1_guess])
    return theta3, l1
    
def velocityCalcs(l2, l3, theta2, theta3, w2, w3_guess=1, v1_guess=1):
    def vector_loop(x):
        w3, v1 = x
        loop1 = -w2*l2*np.sin(theta2)-w3*l3*np.sin(theta3)-v1
        loop2 = w2*l2*np.cos(theta2)+w3*l3*np.cos(theta3)
        return[loop1, loop2]

    try:
        w3, v1 = fsolve(vector_loop, [w3_guess, v1_guess])
        return w3, v1
    except:
        return np.nan, np.nan

def accelerationCalcs(l2, l3, theta2, theta3, w2, w3, a2, a3_guess=1, aL1_guess=1):
    def vector_loop(x):
        a3, aL1 = x
        loop1 = -a2*l2*np.sin(theta2)-a3*l3*np.sin(theta3)-(w2**2)*l2*np.cos(theta2)-(w3**2)*l3*np.cos(theta3)-aL1
        loop2 = a2*l2*np.cos(theta2)+a3*l3*np.cos(theta3)-(w2**2)*l2*np.sin(theta2)-(w3**2)*l3*np.sin(theta3)
        return[loop1, loop2]

    try:
        a3, aL1 = fsolve(vector_loop, [a3_guess, aL1_guess])
        return a3, aL1 
    except:
        return np.nan, np.nan

#settings:
w2_const = 1
a2_const = 0
l2 = 0.105
l3 = 0.172
l4 = 0.027

#position
theta2 = np.linspace((1/12)*np.pi, (1/3)*np.pi, 1000)
theta3=[]
l1=[]

for i in range(len(theta2)):
    temp_theta3, temp_l1 = positionCalcs(l2, l3, l4, theta2[i])
    theta3.append(temp_theta3)
    l1.append(temp_l1)

#velocity
w2=[]
w3=[]
v1=[]
for i in range(len(theta2)):
    w2.append(w2_const)
    temp_w3, temp_v1 = velocityCalcs(l2, l3, theta2[i], theta3[i], w2[i])
    w3.append(temp_w3)
    v1.append(temp_v1)

#acceleration
a2=[]
a3=[]
aL1=[]
for i in range(len(theta2)):
    a2.append(a2_const)
    temp_a3, temp_aL1 = accelerationCalcs(l2, l3, theta2[i], theta3[i], w2[i], w3[i], a2[i])
    a3.append(temp_a3)
    aL1.append(temp_aL1)

# print(theta2)
# print(theta3)
# print(l1)
theta2_deg = [x*(180/np.pi) for x in theta2]

plt.axes()
plt.plot(theta2_deg, aL1,'b-')
plt.xlabel('Theta2 (degrees)')
plt.ylabel('Linear Acceleration of Output (m/s^2)')
plt.title(f'Alpha = {a2_const}, Omega = {w2_const}')
plt.show()