import matplotlib.pyplot as plt
import numpy as np

class PosAnal4Bar:
    def __init__(self, L1, L2, L3, L4, dist, L1_angle=0, dist_theta=0, NumPoints=10000):
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
    
    

# Q2b = PosAnal4Bar(2.22, 0.86, 1.85, 0.86, 1.33)  

# Q2b.graphAngles()
# Q2b.graphP()

Q2d = PosAnal4Bar(2.22, 1.0, 2.06, 2.44, 3.06, L1_angle=(27*np.pi/180), dist_theta = (31*np.pi/180))
Q2d.graphP2()
Q2da = PosAnal4Bar(2.5, 1.0, 2.06, 2.44, 3.06, L1_angle=(27*np.pi/180), dist_theta = (31*np.pi/180))
Q2da.graphP2()
# theQuickNDirty(0.075, 0.17, 0.045, NumPoints=10000, cycles = 4)
