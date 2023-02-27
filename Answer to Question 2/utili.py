import numpy as np
import pandas as pd

'''
Method:  
(1) Define the function of e, s, c, and p.
(2) Define the function of the fourth-order Runge-Kutta method.

f1 = ds/dt, f2 = dc/dt, f3 = de/dt, f4 = dp/dt
*Note: The rates of enzymes are denoted by alpha 1 to 3 due to avoiding notation abuse.
'''

def func_f1(alpha1, alpha2, alpha3, s, e, c):
    return (-alpha1*s*e + alpha2*c)

def func_f2(alpha1, alpha2, alpha3, s, e, c):
    return (alpha1*s*e - alpha2*c -alpha3*c)

def func_f3(alpha1, alpha2, alpha3, s, e, c):
    return (-alpha1*s*e + alpha2*c + alpha3*c)

def func_f4(alpha1, alpha2, alpha3, s, e, c):
    return alpha3 * c

def fn_list(alpha1, alpha2, alpha3, tup):
    s = tup[0]
    c = tup[1]
    e = tup[2]
    p = tup[3]

    return np.array(
        (
            func_f1(alpha1, alpha2, alpha3, s, e, c),
            func_f2(alpha1, alpha2, alpha3, s, e, c),
            func_f3(alpha1, alpha2, alpha3, s, e, c),
            func_f4(alpha1, alpha2, alpha3, s, e, c),
        )
    )

def RK4(h, t_0, s_0, c_0, e_0, p_0, alpha1, alpha2, alpha3):
    args = (t_0, s_0, c_0, e_0, p_0)
    paras = (s_0, c_0, e_0, p_0)

    k1 = fn_list(alpha1, alpha2, alpha3,
        (s_0, c_0, e_0, p_0)
    )  # array of k1 koefficients ( each equation of the system has its separate k1,k2,k3,k4 koefficients)

    k2args = tuple(
        i[0] + h/2 * i[1] / 2 for i in zip(paras, k1)
    )  # array of args for calculating k2
    k2 = fn_list(alpha1, alpha2, alpha3, k2args)  # array of k2

    k3args = tuple(
        i[0] + h/2 * i[1] / 2 for i in zip(paras, k2)
    )  # args for k3
    k3 = fn_list(alpha1, alpha2, alpha3, k3args)  #  array of k3

    k4args = tuple(
        i[0] + i[1] * h for i in zip(paras, k3)
    )  #  args for k4
    k4 = fn_list(alpha1, alpha2, alpha3, k4args)  # array of k4

    new_vals = list(
        args
    )  # values of functions in the next step t = t_0 + h

    new_vals[0] = new_vals[0] + h  # update the t value
    for i in range(1, 5):
        new_vals[i] = new_vals[i] + (k1[i-1] + 2 * k2[i-1] + 2 * k3[i-1] + k4[i-1]) * h / 6  # update function values at the point t = t_0 + h
    return new_vals


def PROCEDURE(verh, h, t_0, s_0, c_0, e_0, p_0, alpha1, alpha2, alpha3):
    steps = (
        int((verh - t_0) / h) + 1
    )  # calculating number of steps required to achieve boundary
    rez = np.zeros((steps, 5))  # results are written inti this array
    for i in range(steps):
        (t_0, s_0, c_0, e_0, p_0) = RK4(h, t_0, s_0, c_0, e_0, p_0, alpha1, alpha2, alpha3)
        rez[i, 0] = t_0
        rez[i, 1] = s_0
        rez[i, 2] = c_0
        rez[i, 3] = e_0
        rez[i, 4] = p_0
    return rez