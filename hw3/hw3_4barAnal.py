import matplotlib.pyplot as plt
import numpy as np

class PosAnal4Bar:
    def __init__(self, L1, L2, L3, L4, L5, L1_angle=0, angularSpeed2=-1, NumPoints=10000):
        self.d = L1
        self.d_angle = L1_angle
        self.a = L2
        self.b = L3
        self.c = L4
        self.L5 = L5
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
        self.theta5=[]
        self.theta4cross=[]
        self.theta3cross=[]
        self.invalidAngles=[]
        self.w2=[]
        self.w3=[]
        self.w4=[]
        self.w5=[]
        self.cXpos=[]
        self.cXvel=[]

        for i in range(self.NumPoints):
            self.theta2.append(i*self.deltaTheta)
            self.w2.append(angularSpeed2)
            self.w3.append(0)
            self.w4.append(0)
            self.cXpos.append(0)
            self.cXvel.append(0)
            self.theta5.append(0)
            self.w5.append(0)

        #calculate values
        self.computeState()
        self.velocityAnalysis()

    def computeState(self):
        for theta in self.theta2:
            A=np.cos(theta)-self.K1-self.K2*np.cos(theta)+self.K3
            B=E=-2*np.sin(theta)
            C=self.K1-(self.K2+1)*np.cos(theta)+self.K3
            D=np.cos(theta)-self.K1+self.K4*np.cos(theta)+self.K5
            F=self.K1+(self.K4-1)*np.cos(theta)+self.K5
            self.updateAngles(A, B, C, D, E, F, theta)

    def updateAngles(self, A, B, C, D, E, F, theta):   
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
            theta4val=2*np.arctan2((-B-np.sqrt(B**2-4*A*C))/(2*A),1)
            self.theta4.append(theta4val)
            theta3val=2*np.arctan2((-E-np.sqrt(E**2-4*D*F))/(2*D),1)
            self.theta3.append(theta3val)
        
            theta4crossval=2*np.arctan2((-B+np.sqrt(B**2-4*A*C))/(2*A),1)
            theta3crossval=2*np.arctan2((-E+np.sqrt(E**2-4*D*F))/(2*D),1)
            self.theta4cross.append(theta4crossval)
            self.theta3cross.append(theta3crossval)

            xpos=self.a*np.cos(theta)
            ypos=self.a*np.sin(theta)
            xposcross=self.a*np.cos(theta)
            yposcross=self.a*np.sin(theta)
            self.xcross.append(xposcross)
            self.ycross.append(yposcross)
            self.x.append(xpos)
            self.y.append(ypos)
            self.xA.append(self.a*np.cos(theta))
            self.yA.append(self.a*np.sin(theta))
            self.xB.append(self.d*np.cos(self.d_angle)+self.c*np.cos(theta4val))
            self.yB.append(self.d*np.sin(self.d_angle)+self.c*np.sin(theta4val))
            self.xB.append(self.d*np.cos(self.d_angle)+self.c*np.cos(theta4crossval))
            self.yB.append(self.d*np.sin(self.d_angle)+self.c*np.sin(theta4crossval))
    
    def velocityAnalysis(self):
        for i in range(len(self.w2)):
            self.w3[i]=(self.a*self.w2[i]/self.b)*(np.sin(self.theta4[i]-self.theta2[i])/np.sin(self.theta3[i]-self.theta4[i]))
            self.w4[i]=(self.a*self.w2[i]/self.c)*(np.sin(self.theta2[i]-self.theta3[i])/np.sin(self.theta4[i]-self.theta3[i]))
            self.theta5[i]=np.arcsin(-self.c*np.sin(self.theta4[i])/self.L5)
            self.cXpos[i]=self.c*np.cos(self.theta4[i])+self.L5*np.cos(self.theta5[i])
            self.w5[i]=(-self.c*self.w4[i]*np.cos(self.theta4[i]))/(self.L5*np.cos(self.theta5[i]))
            self.cXvel[i]=-self.c*self.w4[i]*np.sin(self.theta4[i])-self.L5*self.w5[i]*np.sin(self.theta5[i])


    def graph(self):
        plt.axes().set_aspect('equal')
        plt.plot(self.xA,self.yA,'b-',label='Joint A')
        plt.plot(self.xB,self.yB,'g.',markersize=0.2,label='Joint B')
        plt.legend()
        plt.ylabel('y')
        plt.xlabel('x')
        plt.show()

    def graphFirst(self):
        tempTheta2 = [x*180/np.pi for x in self.theta2]
        plt.plot(tempTheta2,self.cXvel,'b-', label = 'slider Velocity')
        plt.plot(tempTheta2,self.w4,'g-', label = 'angular velocity 4')
        plt.legend()
        plt.ylabel('angularVelocity/linearVelocity')
        plt.xlabel('angle theta2')
        plt.title('Link 4 Angular Velocity and Slider linear Velocity versus Theta 2')
        plt.show()

    def graphSecond(self):
        tempTheta2 = [x*180/np.pi for x in self.theta2]
        plt.plot(tempTheta2,self.cXvel,'b-', label = 'slider Velocity')
        plt.legend()
        plt.ylabel('linearVelocity')
        plt.xlabel('angle theta2')
        plt.title('Slider linear Velocity versus Theta 2')
        plt.show()
    
    def graphThird(self):
        plt.plot(self.cXpos,self.cXvel,'b-', label = 'slider Velocity')
        plt.legend()
        plt.ylabel('linearVelocity')
        plt.xlabel('linearPosition')
        plt.title('Slider linear Velocity versus linear Position')
        plt.show()

    def graphSlider(self):
        tempTheta4 = [x*180/np.pi for x in self.theta4]
        plt.plot(self.w4, self.cXvel,'b-')
        plt.legend()
        plt.ylabel('Velocity of slider')
        plt.xlabel('position of slider')
        plt.show()

    def graphTemp(self):
        tempTheta2 = [x*180/np.pi for x in self.theta2]
        tempTheta4 = [x*180/np.pi for x in self.theta4]
        tempTheta5 = [x*180/np.pi for x in self.theta5]
        plt.plot(tempTheta2, tempTheta4,'*')
        plt.legend()
        plt.xlabel('theta4')
        plt.ylabel('c pos')
        plt.show()

    def percentDeviation(self, startPos, endPos):
        numPoints = 0
        total = 0
        for i in range(len(self.cXvel)):
            if(startPos<self.theta2[i]<endPos):
                numPoints+=1
                total+=self.cXvel[i]
        average = total/numPoints
        totalDeviation = 0
        for i in range(len(self.cXvel)):
            if(startPos<self.theta2[i]<endPos):
                totalDeviation+=abs(self.cXvel[i]-average)
        averageDeviation = totalDeviation/numPoints
        print(f'average: {average}')
        print(f'average deviation: {averageDeviation}')
        percentDeviation=(averageDeviation-average)/average
        return percentDeviation


Q3c = PosAnal4Bar(1.0, 2.17, 2.067, 2.31, L5 = 5.4, L1_angle=(-102*np.pi/180))
# Q3c.graph()
# Q3c.graphVel()
# Q3c.graphSlider()
print("first range: ")
print(Q3c.percentDeviation(240*np.pi/180, 270*np.pi/180))
print("second range: ")
print(Q3c.percentDeviation(190*np.pi/180, 315*np.pi/180))
Q3c.graphFirst()
Q3c.graphSecond()
Q3c.graphThird()
# Q3c.graphTemp()