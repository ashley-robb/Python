'''
ASHLEY ROBB 
STUDENT ID: 20203465
ENPH 213 LAB 4 
'''

import matplotlib.pyplot as plt # used to plot data from file 
import numpy as np # used for importing data 
import math as m # used for math functions such as pi 
import sympy as sm
from sympy.abc import x , a
from sympy import lambdify
import decimal as d

#%% 1. a)
print('\n\n---------------RESULTS FROM QUESTION 1A : ---------------\n\n')


def bisection( func, x_0 , x_1 , kmax = 200 , tol = 1.e-8): # creating a function for bisection method 
    func_0 = func(x_0)
    for k in range(1, kmax): # loop through values 
        x_2 = (x_0 + x_1) / 2
        func_2 = func(x_2)
        
        if func_0 * func_2 < 0: 
            x_1 = x_2
        else:
            x_0 , func_0 = x_2 , func_2
        
        new_x_2 = (x_0 + x_1 )/2
        diff = abs(new_x_2 - x_2)
        
        if abs(diff / new_x_2) < tol:
            break
    else:
        new_x_2 = None
    return new_x_2

function1a = 1/(x-3) # defining the function 
f1a = lambdify([x], function1a) # lambdifying the function 
root1a = bisection(f1a , 0 , 5) 

print('There is a root at: ' , root1a , '\n') # printing the function 

print('Justification of answer: The outputted answer is not correct because the')
print('given function 1 / (x-3) does not have any real roots, when x = 3, the')
print('solution is not a number (cannot divide by zero')

#%% 1. b)
print('\n\n---------------RESULTS FROM QUESTION 1B : ---------------\n\n')

function1b = sm.exp( x - sm.sqrt(x) ) - x # defining the given function symbolically 
f1b = lambdify([x], function1b) # lambdifying the function 

root1b_bisec = bisection(f1b , 0 , 1.5) 

print('Using bisection, we can determine that there is a root at: ' , root1b_bisec , '\n')


def diff(func, x , tol): # derivative function 
    return (func(x + tol / 2) - func(x - tol/2)) / tol
    

def newton(func , x_0 , kmax = 200 , tol = 1.e-8): # loop for newtons function 
    x = np.zeros(kmax + 1) # array of values 
    x[0] = x_0 # x initial 
    for i in range(0,kmax): # loop thorugh values 
        j = i + 1
        derive = diff(func , x[i] , tol) # derivative 
        x[j] = x[i] - (func(x[i])) / derive # new x 
        if abs(func(x[j])) > tol: # tolerance test 
            plt.plot([x[i] , x[i] ], [func(x[i]) , 0] ,'k:') # plotting results   
            plt.plot([x[j], x[i]] , [0 , func(x[i])] , 'r--')
        else: 
            break 
    return x[j] # returning desired value 

root1b_newt = newton(f1b , 0.01)
            
print("Using Newton's method, we determine that there is a root at : " , root1b_newt)

x_vals = np.linspace(0, 1.1, 200) # x values for plotting 
y_vals = f1b(x_vals) # y values for plotting 

plt.plot(x_vals , y_vals ,'b-') # plotting the function 
plt.axhline(0, color = 'black' , linewidth = 0.75) # adding the axis line at y = 0
plt.xlabel('x' , size = 20) ; plt.ylabel('f(x)' , size = 20) # labeling the axis'
plt.grid() # adding grid background 
plt.annotate('$x^{(0)}$' , [-0.05 , -0.09] , size = 15) # adding x_0 label 
plt.annotate('$x^{(1)}$' , [0.15 , -0.09] , size = 15) # adding x_1 label 
plt.annotate('$x^{(2)}$' , [0.68 , -0.09] , size = 15) # x_2 label 
plt.annotate('$x^{(3)}$' , [0.93 , -0.09] , size = 15) # x_3 label


#%% 1. c)
print('\n\n---------------RESULTS FROM QUESTION 1C : ---------------\n\n')


def newton_2nd_root(func , x , a_val ,kmax = 200 , tol = 1.e-8):  # defining a function for the second root 
    for i in range(1,kmax): # loop 
        u = func(x)/(x-a_val) # defining u 
        u_prime = (diff(func , x , tol) * (x-a_val) - func(x)) * ((x - a_val)**(-2)) # new u 
        if abs(u_prime) > tol: # tolerance test
            x = x - u/u_prime # new x 
        else: 
            break 
    return x # desired output 

# printing outputs for necessary test points 

print('2nd root using starting point of 2:' , newton_2nd_root(f1b, 2, 1))
print('2nd root using starting point of 0.5:' , newton_2nd_root(f1b, 2, 1))
print('2nd root using starting point of 4:' , newton_2nd_root(f1b, 2, 1))
print('2nd root using starting point of 0.1:' , newton_2nd_root(f1b, 2, 1))














