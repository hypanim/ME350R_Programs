from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def four_bar(link_lengths:List[float], input_angle:float, theta_1:float = 0.0, angle_p:float = 0.0, dist_p:float = None, cross_solution:bool = False) -> Tuple[float]:
    #four bar linkage equations
    if dist_p is None:
        dist_p = link_lengths[2]

    k1 = link_lengths[0]/link_lengths[1]
    k2 = link_lengths[0]/link_lengths[3]
    k3 = (link_lengths[1]**2 - link_lengths[2]**2 + link_lengths[3]**2 + link_lengths[0]**2) / (2*link_lengths[1]*link_lengths[3])
    k4 = link_lengths[0]/link_lengths[2]
    k5 = (link_lengths[3]**2 - link_lengths[0]**2 - link_lengths[1]**2 - link_lengths[2]**2) / (2*link_lengths[1]*link_lengths[2])

    A = np.cos(input_angle) - k1 - k2*np.cos(input_angle) + k3
    B = -2*np.sin(input_angle)
    C = k1 - (k2+1)*np.cos(input_angle) + k3
    D = np.cos(input_angle) - k1 + k4*np.cos(input_angle) + k5
    E = -2*np.sin(input_angle)
    F = k1 + (k4-1)*np.cos(input_angle) + k5

    flip = -1 if cross_solution else 1
    theta_4 = 2*np.arctan((-B + flip*np.sqrt(B**2 - 4*A*C))/(2*A))
    theta_3 = 2*np.arctan((-E + flip*np.sqrt(E**2 - 4*D*F))/(2*D))

    input_angle += theta_1
    theta_3 += theta_1
    theta_4 += theta_1

    p_x = link_lengths[1]*np.cos(input_angle) + dist_p*np.cos(theta_3+angle_p)
    p_y = link_lengths[1]*np.sin(input_angle) + dist_p*np.sin(theta_3+angle_p)

    return theta_3, theta_4, p_x, p_y

def four_bar_fs(link_lengths:List[float], input_angle:float, theta_1:float = 0.0, angle_p:float = 0.0, dist_p:float = None, theta_3init = 0, theta_4init = 0) -> Tuple[float]:
#Equations for a four-bar linkage using f-solve
    if dist_p is None:
        dist_p = link_lengths[2]

    def vector_loop(x):
        #vector loop
        theta_3, theta_4 = x
        # r1 * cos(theta_1) + r2 * cos(theta_2) ...
        loop1 = link_lengths[1]*np.cos(input_angle) + link_lengths[2]*np.cos(theta_3) - link_lengths[3]*np.cos(theta_4) - link_lengths[0]
        loop2 = link_lengths[1]*np.sin(input_angle) + link_lengths[2]*np.sin(theta_3) - link_lengths[3]*np.sin(theta_4)
        return [loop1, loop2]

    theta_3, theta_4 = fsolve(vector_loop, [theta_3init, theta_4init])

    input_angle += theta_1
    theta_3 += theta_1
    theta_4 += theta_1

    p_x = link_lengths[1]*np.cos(input_angle) + dist_p*np.cos(theta_3+angle_p)
    p_y = link_lengths[1]*np.sin(input_angle) + dist_p*np.sin(theta_3+angle_p)

    return theta_3, theta_4, p_x, p_y

def main():
#main function
    ts = np.linspace(-np.pi/4, (1/4)*np.pi, 50)
    xs = []
    ys = []
    theta3 = []
    theta4 = []
    link_lengths = [2.22, 0.86, 1.85, 0.86]

    #used for getting initial guess
    theta3_guess, theta4_guess, p_x, P_y = four_bar_fs(link_lengths, 0, theta_1 = np.deg2rad(0), angle_p = 0, dist_p = 1.33, theta_3init = 0.84, theta_4init = 2)
    theta3.append(theta3_guess)
    theta4.append(theta4_guess)
    print(f'theta 3 is {theta3}, theta 4 is {theta4}')

    for input_angle in ts:
        theta3_guess, theta4_guess, p_x, p_y = four_bar_fs(link_lengths, input_angle, theta_1 = np.deg2rad(0), angle_p = 0, dist_p = 1.33, theta_3init=theta3[-1], theta_4init=theta4[-1])
        xs.append(p_x)
        ys.append(p_y)
        theta3.append(theta3_guess)
        theta4.append(theta4_guess)
        print(f'for input angle of {input_angle*(180/np.pi)}, theta3 = {theta3[-1]*(180/np.pi)}, theta4 = {theta4[-1]*(180/np.pi)}')
    plt.plot(xs, ys, '.')
    plt.show()

if __name__ == '__main__':
    main()