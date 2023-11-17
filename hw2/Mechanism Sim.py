
import matplotlib.pyplot as plt
import numpy as np
import math as m

NumPoints=10000
deltaTheta=2*np.pi/NumPoints

# System Geometry Parameters
L1=d=2.22
L2=a=0.86
L3=b=1.85
L4=c=0.86
dist=1.33

#Initialise Position and Angle vectors
x=[]
y=[]
xcross=[]
ycross=[]
xA=[]
yA=[]
xB=[]
yB=[]
theta2=[]
theta4=[]
theta3=[]
theta4cross=[]
theta3cross=[]
invalidAngles=[]

# Compute Constant Coefficients
K1=d/a 
K2=d/c 
K3=(a**2-b**2+c**2+d**2)/(2*a*c)
K4=d/b
K5=(c**2-d**2-a**2-b**2)/(2*a*b)

# Generate input angle array
for i in range(NumPoints):
    theta2.append(i*deltaTheta)

# Compute Linkage State for all angles
for theta in theta2:
    A=np.cos(theta)-K1-K2*np.cos(theta)+K3
    B=E=-2*np.sin(theta)
    C=K1-(K2+1)*np.cos(theta)+K3
    D=np.cos(theta)-K1+K4*np.cos(theta)+K5
    F=K1+(K4-1)*np.cos(theta)+K5

    if B**2-4*A*C < 0 or E**2-4*D*F < 0:
        theta4.append(np.nan)
        theta3.append(np.nan)
        xcross.append(np.nan)
        ycross.append(np.nan)
        x.append(np.nan)
        y.append(np.nan)
        xA.append(np.nan)
        yA.append(np.nan)
        xB.append(np.nan)
        yB.append(np.nan)
        invalidAngles.append(theta)
    else:
        theta4val=2*np.arctan2((-B-np.sqrt(B**2-4*A*C))/(2*A),1)
        theta4.append(theta4val)
        theta3val=2*np.arctan2((-E-np.sqrt(E**2-4*D*F))/(2*D),1)
        theta3.append(theta3val)
    
        theta4crossval=2*np.arctan2((-B+np.sqrt(B**2-4*A*C))/(2*A),1)
        theta3crossval=2*np.arctan2((-E+np.sqrt(E**2-4*D*F))/(2*D),1)
        theta4cross.append(theta4crossval)
        theta3cross.append(theta3crossval)

        xpos=a*np.cos(theta)+dist*np.cos(theta3val)
        ypos=a*np.sin(theta)+dist*np.sin(theta3val)
        xposcross=a*np.cos(theta)+dist*np.cos(theta3crossval)
        yposcross=a*np.sin(theta)+dist*np.sin(theta3crossval)
        xcross.append(xposcross)
        ycross.append(yposcross)
        x.append(xpos)
        y.append(ypos)
        xA.append(a*np.cos(theta))
        yA.append(a*np.sin(theta))
        xB.append(d+c*np.cos(theta4val))
        yB.append(c*np.sin(theta4val))
        xB.append(d+c*np.cos(theta4crossval))
        yB.append(c*np.sin(theta4crossval))

plt.axes().set_aspect('equal')
plt.plot(x,y,'r-',label='Point P')
plt.plot(xA,yA,'b-',label='Joint A')
plt.plot(xB,yB,'g.',markersize=0.2,label='Joint B')
plt.plot(0,0,'bo')
plt.plot(d,0,'go')
plt.plot(xcross,ycross,'r--',label='Point P Cross Angle')
plt.legend()
plt.ylabel('y')
plt.xlabel('x')
plt.show()


"""
plt.plot([1, 2, 3, 4])
plt.ylabel('some numbers')
plt.show()
"""