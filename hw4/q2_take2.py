import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

class PosAnal4Bar:
    def __init__(self, L1, L2, L3, L4, dist, L1_angle=0, dist_theta=0, NumPoints=1000):
        self.d = L1
        self.d_angle = L1_angle
        self.a = L2
        self.b = L3
        self.c = L4
        self.dist = dist
        self.dist_theta = dist_theta
        self.NumPoints = NumPoints

        self.deltaTheta=2*np.pi/NumPoints

        self.K1=self.d/self.a 
        self.K2=self.d/self.c 
        self.K3=(self.a**2-self.b**2+self.c**2+self.d**2)/(2*self.a*self.c)
        self.K4=self.d/self.b
        self.K5=(self.c**2-self.d**2-self.a**2-self.b**2)/(2*self.a*self.b)

        #initialize lists
        self.x=[]
        self.y=[]
        self.xcross=[]
        self.ycross=[]
        self.xA=[]
        self.yA=[]
        self.xB=[]
        self.yB=[]
        self.theta2=[]
        self.theta4=[]
        self.theta3=[]
        self.theta4cross=[]
        self.theta3cross=[]
        self.invalidAngles=[]

        for i in range(self.NumPoints):
            self.theta2.append(i*self.deltaTheta)

        #calculate values
        self.computeState()

    def computeState(self):
        for theta in self.theta2:
            A=np.cos(theta)-self.K1-self.K2*np.cos(theta)+self.K3
            B=E=-2*np.sin(theta)
            C=self.K1-(self.K2+1)*np.cos(theta)+self.K3
            D=np.cos(theta)-self.K1+self.K4*np.cos(theta)+self.K5
            F=self.K1+(self.K4-1)*np.cos(theta)+self.K5
            self.updateAngles(A, B, C, D, E, F, theta, self.dist)

    def updateAngles(self, A, B, C, D, E, F, theta, dist):   
        if B**2-4*A*C < 0 or E**2-4*D*F < 0:
            self.theta4.append(np.nan)
            self.theta4cross.append(np.nan)
            self.theta3.append(np.nan)
            self.theta3cross.append(np.nan)
            self.xcross.append(np.nan)
            self.ycross.append(np.nan)
            self.x.append(np.nan)
            self.y.append(np.nan)
            self.xA.append(np.nan)
            self.yA.append(np.nan)
            self.xB.append(np.nan)
            self.yB.append(np.nan)
            self.invalidAngles.append(theta)
        else:
            theta4val=2*np.arctan2((-B-np.sqrt(B**2-4*A*C))/(2*A),1)+self.d_angle
            self.theta4.append(theta4val)
            theta3val=2*np.arctan2((-E-np.sqrt(E**2-4*D*F))/(2*D),1)+self.d_angle
            self.theta3.append(theta3val)
        
            theta4crossval=2*np.arctan2((-B+np.sqrt(B**2-4*A*C))/(2*A),1)+self.d_angle
            theta3crossval=2*np.arctan2((-E+np.sqrt(E**2-4*D*F))/(2*D),1)+self.d_angle
            self.theta4cross.append(theta4crossval)
            self.theta3cross.append(theta3crossval)

            xpos=self.a*np.cos(theta+self.d_angle)+dist*np.cos(theta3val+self.dist_theta)
            ypos=self.a*np.sin(theta+self.d_angle)+dist*np.sin(theta3val+self.dist_theta)
            xposcross=self.a*np.cos(theta+self.d_angle)+dist*np.cos(theta3crossval+self.dist_theta)
            yposcross=self.a*np.sin(theta+self.d_angle)+dist*np.sin(theta3crossval+self.dist_theta)
            self.xcross.append(xposcross)
            self.ycross.append(yposcross)
            self.x.append(xpos)
            self.y.append(ypos)
            self.xA.append(self.a*np.cos(theta+self.d_angle))
            self.yA.append(self.a*np.sin(theta+self.d_angle))
            self.xB.append(self.d*np.cos(self.d_angle)+self.c*np.cos(theta4val))
            self.yB.append(self.d*np.sin(self.d_angle)+self.c*np.sin(theta4val))
            self.xB.append(self.d*np.cos(self.d_angle)+self.c*np.cos(theta4crossval))
            self.yB.append(self.d*np.sin(self.d_angle)+self.c*np.sin(theta4crossval))
    
    def graph(self):
        plt.axes().set_aspect('equal')
        plt.plot(self.x,self.y,'r-',label='Point P')
        plt.plot(self.xA,self.yA,'b-',label='Joint A')
        plt.plot(self.xB,self.yB,'g.',markersize=0.2,label='Joint B')
        plt.plot(0,0,'bo')
        plt.plot(self.d,0,'go')
        plt.plot(self.xcross,self.ycross,'r--',label='Point P Cross Angle')
        plt.legend()
        plt.ylabel('y')
        plt.xlabel('x')
        plt.show()

    def graphAngles(self):
        theta2_deg = [x*(180/np.pi) for x in self.theta2]
        theta3_deg = [x*(180/np.pi) for x in self.theta3]
        theta3Cross_deg = [x*(180/np.pi) for x in self.theta3cross]

        plt.axes().set_aspect('equal')
        plt.plot(theta2_deg,theta3_deg,'b-',label='Theta 3 angle')
        plt.plot(theta2_deg,theta3Cross_deg, 'g--',label='Theta 3 inverted angle')
        plt.legend()
        plt.title("Graph of Joint Angles")
        plt.xlabel('Angle (degrees)')
        plt.ylabel('Angle (degrees)')
        plt.show()

    def graphP(self):
        plt.axes().set_aspect('equal')
        plt.plot(self.x,self.y,'b-',label='Point P')
        plt.plot(self.xcross,self.ycross,'g--',label='Point P inverted')
        plt.legend()
        plt.title("Graph of Point P")
        plt.xlabel('X (in)')
        plt.ylabel('Y (in)')
        plt.show()
        
    def graphP2(self):
        plt.axes().set_aspect('equal')
        plt.plot(self.xcross,self.ycross,'b-',label='Point P')
        plt.legend()
        plt.title("Graph of Point P")
        plt.xlabel('X (in)')
        plt.ylabel('Y (in)')
        plt.show()

def theQuickNDirty(L2, L3, h, NumPoints=10000, cycles=1):
    deltaTheta=2*np.pi*cycles/NumPoints
    theta2=[]
    thetaB=[]
    Bx=[]
    for i in range(NumPoints):
        theta2.append(i*deltaTheta)

    for i in range(NumPoints):
        thetaB.append(np.arcsin((L2*np.sin(theta2[i])-h)/L3))
    
    for i in range(NumPoints):
        Bx.append(-L2*np.cos(theta2[i])+L3*np.cos(thetaB[i]))

    theta2_deg = [x*(180/np.pi) for x in theta2]
    plt.axes()
    plt.plot(theta2_deg, Bx,'b-',label='Bx')
    plt.legend()
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Bx (m)')
    plt.show()
    

Q2b = PosAnal4Bar(2.22, 0.86, 1.85, 0.86, 1.33)  
# Q2b.graphAngles()
# Q2b.graphP2()

theta2 = Q2b.theta2.copy()
theta3 = Q2b.theta3cross.copy()
theta4 = Q2b.theta4cross.copy()

#cleaning data 
while np.nan in theta3:
    removalIndex = theta3.index(np.nan)
    theta2.pop(removalIndex)
    theta3.pop(removalIndex)
    theta4.pop(removalIndex)

#velocity analysis:

def velocityCalcs(theta2, theta3, theta4, w2, l2, l3, l4, w3_guess=0, w4_guess=0):
    def vector_loop(x):
        w3, w4 = x
        loop1 = -l2*w2*np.sin(theta2)-l3*w3*np.sin(theta3)+l4*w4*np.sin(theta4)
        loop2 = l2*w2*np.cos(theta2) + l3*w3*np.cos(theta3) - l4*w4*np.cos(theta4)
        return [loop1, loop2]

    try:
        w3, w4 = fsolve(vector_loop, [w3_guess, w4_guess])
        return w3, w4
    except:
        return np.nan, np.nan

w2 = []
w3 = []
w4 = []

for i in range(len(theta2)):
    w2.append(1)
    temp_w3, temp_w4 = velocityCalcs(theta2[i], theta3[i], theta4[i], w2[i], l2 = 0.86, l3 = 1.85, l4 = 0.86)
    w3.append(temp_w3)
    w4.append(temp_w4)

# plt.axes()
# plt.plot(theta2, w2, '.')
# plt.plot(theta2, w3, '.')
# plt.plot(theta2, w4, '.')
# plt.show()

#acceleration analysis
def accelerationCalcs(theta2, theta3, theta4, w2, w3, w4, l2, l3, l4, a2, a3_guess=0, a4_guess=0):
    def vector_loop(x):
        a3, a4 = x
        loop1 = -a4*l4*np.cos(theta4)+a3*l3*np.cos(theta3)+a2*l2*np.cos(theta2)+(w4**2)*l4*np.sin(theta4)-(w3**2)*l3*np.sin(theta3)-(w2**2)*l2*np.sin(theta2)
        loop2 = (w4**2)*l4*np.cos(theta4)-(w3**2)*l3*np.cos(theta3)-(w2**2)*l2*np.cos(theta2)+a4*l4*np.sin(theta4)-a3*l3*np.sin(theta3)-a2*l2*np.sin(theta2)
        return [loop1, loop2]

    try:
        a3, a4 = fsolve(vector_loop, [a3_guess, a4_guess])
        return a3, a4
    except:
        return np.nan, np.nan

a2 = []
a3 = []
a4 = []

for i in range(len(theta2)):
    a2.append(1)
    temp_a3, temp_a4 = accelerationCalcs(theta2[i], theta3[i], theta4[i], w2[i], w3[i], w4[i], 0.86, 1.85, 0.86, a2[i])
    a3.append(temp_a3)
    a4.append(temp_a4)

# theta2_deg = [x*(180/np.pi) for x in theta2]
# plt.axes()
# plt.plot(theta2_deg, a3, '.')
# plt.xlabel('Theta2 (degrees)')
# plt.ylabel('Angular acceleration of the coupler (rad/s^2)')
# plt.title('Acceleration of the coupler for Crossed')
# plt.ylim([-5, 5])
# plt.show()

#part C
theta2_uc = Q2b.theta2.copy()
theta2_c = Q2b.theta2.copy()
theta3_uc = Q2b.theta3.copy()
theta3_c = Q2b.theta3cross.copy()
theta4_uc = Q2b.theta4.copy()
theta4_c = Q2b.theta4cross.copy()

#cleaning data 
while np.nan in theta3_uc:
    removalIndex = theta3_uc.index(np.nan)
    theta2_uc.pop(removalIndex)
    theta3_uc.pop(removalIndex)
    theta4_uc.pop(removalIndex)

#cleaning data 
while np.nan in theta3_c:
    removalIndex = theta3_c.index(np.nan)
    theta2_c.pop(removalIndex)
    theta3_c.pop(removalIndex)
    theta4_c.pop(removalIndex)

#uncrossed first
w2_uc = []
w3_uc = []
w4_uc = []

for i in range(len(theta2_uc)):
    w2_uc.append(1)
    temp_w3_uc, temp_w4_uc = velocityCalcs(theta2_uc[i], theta3_uc[i], theta4_uc[i], w2_uc[i], l2 = 0.86, l3 = 1.85, l4 = 0.86)
    w3_uc.append(temp_w3_uc)
    w4_uc.append(temp_w4_uc)

a2_uc = []
a3_uc = []
a4_uc = []

for i in range(len(theta2_uc)):
    a2_uc.append(1)
    temp_a3_uc, temp_a4_uc = accelerationCalcs(theta2_uc[i], theta3_uc[i], theta4_uc[i], w2_uc[i], w3_uc[i], w4_uc[i], 0.86, 1.85, 0.86, a2_uc[i])
    a3_uc.append(temp_a3_uc)
    a4_uc.append(temp_a4_uc)

#now crossed
w2_c = []
w3_c = []
w4_c = []

for i in range(len(theta2_c)):
    w2_c.append(1)
    temp_w3_c, temp_w4_c = velocityCalcs(theta2_c[i], theta3_c[i], theta4_c[i], w2_c[i], l2 = 0.86, l3 = 1.85, l4 = 0.86)
    w3_c.append(temp_w3_c)
    w4_c.append(temp_w4_c)

a2_c = []
a3_c = []
a4_c = []

for i in range(len(theta2_c)):
    a2_c.append(1)
    temp_a3_c, temp_a4_c = accelerationCalcs(theta2_c[i], theta3_c[i], theta4_c[i], w2_c[i], w3_c[i], w4_c[i], 0.86, 1.85, 0.86, a2_c[i])
    a3_c.append(temp_a3_c)
    a4_c.append(temp_a4_c)

# plt.axes()
# plt.plot(theta2_uc, a2_uc, '.')
# plt.plot(theta2_uc, a3_uc, '.')
# plt.plot(theta2_uc, a4_uc, '.')
# plt.plot(theta2_c, a2_c, '.')
# plt.plot(theta2_c, a3_c, '.')
# plt.plot(theta2_c, a4_c, '.')
# plt.ylim([-5, 5])
# plt.show()

theta2_deg = [x*(180/np.pi) for x in theta2]
acceleration_p=[]
p_length = 1.33
for i in range(len(theta2)):
    tangential_i = -a3[i]*p_length*np.sin(theta3[i])
    tangential_j =a3[i]*p_length*np.cos(theta3[i])
    radial_i = -(w3[i]**2)*p_length*np.cos(theta3[i])
    radial_j = -(w3[i]**2)*p_length*np.sin(theta3[i])
    i_total = tangential_i+radial_i
    j_total = tangential_j+radial_j
    total = np.sqrt(i_total**2+j_total**2)
    acceleration_p.append(total)

plt.axes()
plt.plot(theta2_deg, acceleration_p, '.')
plt.xlabel('Theta2 (degrees)')
plt.ylabel('Angular acceleration of the coupler (rad/s^2)')
plt.title('Instaneous Magnitude of Acceleration at P (m/s^2)')
plt.ylim([-5, 5])
plt.show()