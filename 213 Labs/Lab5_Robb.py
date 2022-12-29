"""
213 Lab 4 
Ashley Robb 
03/01/2022
"""

import matplotlib.pyplot as plt # used to plot data from file 
import numpy as np # used for importing data 
import math as m # used for math functions such as pi 
import sympy as sm
import scipy as sc
from sympy import lambdify
from sympy.abc import x 

#%% Q1 PART A

print('----------RESULTS FOR QUESTION 1 PART A: ----------')

def gen_data(n , func , spacing): #creating a function to generate cheb or equidistant data 
    xs = np.zeros(n)
    function = lambdify([x] , func)
    if (spacing == 'chebyshev'): # cheb section 
        for i in range(0,n):
            xs[i] = -m.cos((i*m.pi) / (n-1))
    elif (spacing == 'equidistant'): # equid section
        xs = np.linspace(-1 , 1 , n)
    else:
        return('ERROR: The spacing type entered is not recognized. The spacing input should be chebyshev or equidistant') #error message for wrong inputs
    f = function(xs)
    return(xs , f)

function_runge = 1/(1+25*(x**2)) # defining the runge function 

# setting x and y values used throughout 
xs_cheb = gen_data(15 , function_runge, 'chebyshev')[0]
ys_cheb = gen_data(15 , function_runge, 'chebyshev')[1]

xs_equid = gen_data(15 , function_runge, 'equidistant')[0]
ys_equid = gen_data(15 , function_runge, 'equidistant')[1]

def monomial(x_vals , y_vals , n_points): # defining a function for monomial
    x_matrix = np.vander(x_vals, increasing = True)
    c = sc.linalg.solve(x_matrix , y_vals)
    x_new = np.vander(np.linspace(-1 , 1 , n_points), len(x_vals) , increasing = True)
    y_new = x_new@c
    return(x_new[:,1] , y_new)

def lagrange(x_vals , y_vals , n_points): # defining a function for lagrange 
    x_new = np.linspace(-1, 1 , n_points)
    y_new = np.zeros(n_points)
    for i in range(0 , len(y_vals)):
        p = 1
        for j in range(0 , len(x_vals)):
            if j != i:
                p*= (x_new - x_vals[j]) / (x_vals[i] -  x_vals[j])
        y_new+= p*y_vals[i]
    return x_new , y_new

y_mc = monomial(xs_cheb, ys_cheb , 100)[1] # y monomial cheb
x_mc = monomial(xs_cheb , ys_cheb , 100)[0] # x monomial cheb

y_me = monomial(xs_equid , ys_equid , 100)[1] # y monomial equid
x_me = monomial(xs_equid , ys_equid , 100)[0] # x monomial equid


y_lc = lagrange(xs_cheb, ys_cheb , 100)[1] # y lagrange cheb
x_lc = lagrange(xs_cheb , ys_cheb , 100)[0] # x lagrange cheb 

y_le = lagrange(xs_equid , ys_equid , 100)[1] # y lagrange equid 
x_le = lagrange(xs_equid , ys_equid , 100)[0] # x lagrange equid

# cheb plot
plt.subplot(1,2,1)
plt.plot(xs_cheb , ys_cheb, 'mo' , label = 'points')
plt.plot(x_lc , y_lc , 'r' , label = 'Lagrange')
plt.plot(x_mc , y_mc, 'g--' , label = 'Monomial')
plt.legend(loc = 'best')
plt.ylabel('f')
plt.xlabel('x')


# equid plot
plt.subplot(1,2,2)
plt.plot(xs_equid , ys_equid , 'mo' , label = 'points')
plt.plot(x_le , y_le , 'r' , label = 'Lagrange')
plt.plot(x_me , y_me, 'g--' , label = 'Monomial')
plt.legend(loc = 'best')
plt.xlabel('x')
plt.show()

print('for question 1a results view plots')

#%% Q1 PART B

print('----------RESULTS FOR QUESTION 1 PART B: ----------')

xs_cheb91 = gen_data(91 , function_runge, 'chebyshev')[0] # 91 cheb data
ys_cheb91 = gen_data(91 , function_runge, 'chebyshev')[1]

y_m91 = monomial(xs_cheb91, ys_cheb91 , 5000)[1] # 91 cheb monomial 
x_m91 = monomial(xs_cheb91 , ys_cheb91 , 5000)[0]
y_l91 = lagrange(xs_cheb91, ys_cheb91 , 5000)[1] # 91 cheb lagrange
x_l91 = lagrange(xs_cheb91 , ys_cheb91 , 5000)[0]


xs_cheb101 = gen_data(101 , function_runge, 'chebyshev')[0] # 101 cheb data 
ys_cheb101 = gen_data(101 , function_runge, 'chebyshev')[1]

y_m101 = monomial(xs_cheb101, ys_cheb101 , 5000)[1] # 101 cheb monomial 
x_m101 = monomial(xs_cheb101 , ys_cheb101 , 5000)[0]

y_l101 = lagrange(xs_cheb, ys_cheb , 5000)[1] # 101 cheb lagrange 
x_l101 = lagrange(xs_cheb , ys_cheb , 5000)[0]

# cheb plot 
plt.subplot(1,2,1)
plt.plot(xs_cheb , ys_cheb, 'mo' , label = 'points')
plt.plot(x_l91 , y_l91 , 'r' , label = 'Lagrange')
plt.plot(x_m91 , y_m91, 'g--' , label = 'Monomial')
plt.legend(loc = 'best')
plt.title('n = 91')
plt.ylabel('f')
plt.xlabel('x')


# equid plot
plt.subplot(1,2,2)
plt.plot(xs_equid , ys_equid , 'mo' , label = 'points')
plt.plot(x_l101 , y_l101 , 'r' , label = 'Lagrange')
plt.plot(x_m101 , y_m101, 'g--' , label = 'Monomial')
plt.legend(loc = 'best')
plt.title('n = 101')
plt.xlabel('x')
plt.show()

print('\n\nFINDINGS: I noticed that the lagrange and monomial methods have some differences for the n = 101 case, I also notice that python outputs an error about accuracy.')

print('for question 1b results view plots')

#%% Q2

print('----------RESULTS FOR QUESTION 2: ----------')

def cubic(x1 , x_vals , y_vals): # defining the cubic spline function 
    n = len(x_vals) 
    tridi = np.zeros((n-2 , n-2)) # creating an array for the tridiagonal 
    np.fill_diagonal(tridi , 2*(x_vals[2:]- x_vals[:-2])) # filling the diagonal rows 
    np.fill_diagonal(tridi[:,1:], x_vals[2:-1]- x_vals[1:-2])
    np.fill_diagonal(tridi[1:,:] , x_vals[2:-1] - x_vals[1:-2])
    c = np.zeros(n)
    b = 6*( ((y_vals[2:]- y_vals[1:-1])/ (x_vals[2:] - x_vals[1:-1]) - (y_vals[1:-1] - y_vals[:-2]) / (x_vals[1:-1] - x_vals[:-2])) )
    c[1:-1] = np.linalg.solve(tridi , b)
    above = np.argmax(x_vals > x1) # determining x points above inputted x values 
    x_above = x_vals[above]
    x_above1 = x_vals[above-1]
    y_above = y_vals[above]
    y_above1 = y_vals[above-1]
    c_above = c[above]
    c_above1 = c[above-1]
    curve = y_above1 * ((x_above - x1) / (x_above - x_above1)) + y_above*((x1 - x_above1) / (x_above - x_above1)) 
    curve = curve - (c_above1/6) * ((x_above - x_above1) * (x_above - x_above1) - (x_above -x1)**3 / (x_above - x_above1))
    curve = curve - (c_above/6) * ((x1 - x_above1)* (x_above - x_above1) - (x1 - x_above1)**3 / (x_above - x_above1))
    return curve


xs_cheb7 = gen_data(7 , function_runge, 'chebyshev')[0] # x 7 cheb 
ys_cheb7 = gen_data(7 , function_runge, 'chebyshev')[1] # y7 cheb

xs_equid7 = gen_data(7 , function_runge, 'equidistant')[0] # x 7 equid
ys_equid7 = gen_data(7 , function_runge, 'equidistant')[1] # y 7 equid 

# data for plots 
xs = np.linspace(-1 , 1 , 500)
y_cheb15 = np.zeros(500)
y_cheb7 = np.zeros(500)
y_equid15 = np.zeros(500)
y_equid7 = np.zeros(500)

for j in range(0,500):
    y_cheb15[j] = cubic(xs[j], xs_cheb, ys_cheb)
    y_cheb7[j] = cubic(xs[j] , xs_cheb7 , ys_cheb7)
    y_equid15[j] = cubic(xs[j] , xs_equid , ys_equid)
    y_equid7[j] = cubic(xs[j] , xs_equid7 , ys_equid7)

# cheb plot
plt.subplot(1,2,1)
plt.plot(xs_cheb , ys_cheb, 'mo' , label = 'points')
plt.plot(xs , y_cheb15 , 'r' , label = 'n=15')
plt.plot(xs , y_cheb7, 'g--' , label = 'n=7')
plt.legend(loc = 'best')
plt.ylabel('f')
plt.xlabel('x')

#equid plot
plt.subplot(1,2,2)
plt.plot(xs_equid , ys_equid , 'mo' , label = 'points')
plt.plot(xs , y_equid15 , 'r' , label = 'n=15')
plt.plot(xs , y_equid7, 'g--' , label = 'n=7')
plt.legend(loc = 'best')
plt.xlabel('x')

plt.show()

print('question 2 results graphed above')

#%% Q3

print('----------RESULTS FOR QUESTION 3: ----------')


def func_q3(x): # defining function for question 3
    y = np.exp(np.sin(2*x))
    return y

def ds(n): # data creating function 
    x_vals = np.zeros(n)
    for i in range(0,n):
        x_vals[i] = (2*np.pi * i) / n
    return x_vals

x11 = ds(11)
y11 = func_q3(x11)

x51 = ds(51)
y51 = func_q3(x51)

def trig(x_vals , y_vals): # function for trig method 
    j = int((len(x_vals) - 1) /2)
    n = int(len(x_vals))
    a = []
    for i in range(0 , j+1): 
        b = 0 
        for k in range (n): #  loop through all vals x and y
            s = y_vals[k] * np.cos(i * x_vals[k])
            b += s
        a.append((1/j) * b)
    c = []
    for i in range(j): # loop through (n-1) / 2 values 
        d = 0
        for k in range(n): # loop through all x and y vals 
            s = y_vals[k] * np.sin(i * x_vals[k])
            d += s
        c.append((1/j) * d)
    c.append(0)
    p = 0
    x = sm.symbols('x')
    for i in range(1 , j): # loop from 1 to (n-1) /2
        s = a[i] * sm.cos(i *x) + c[i] * sm.sin(i *x)
        p += s
    q = (1/2) * a[0] +p
    q = sm.lambdify(x,q, 'numpy')
    return q # return output 

q3_11 = trig(x11 , y11) # running vals through trig func       
x500 = ds(500) # new data 
q3_51 = trig(x51 , y51) # new data through  trig function 
y_500 = q3_51(x500) # new data trough given function 

# n = 11 plot (left)
plt.subplot(1 , 2 , 1)
plt.plot(x11 , y11, 'mo' , label = 'points')
plt.plot(x500 , y_500 , 'r' , label = 'n = 11 (trig |')
plt.xlabel('x')
plt.ylabel('f')
plt.legend(loc = 'best')

# n = 51 plot(right)
plt.subplot(1 , 2, 2)     
plt.plot(x51 , y51, 'mo' , label = 'points')
plt.plot(x500 , y_500 , 'r' , label = 'n = 51 (trig |)')
plt.xlabel('x')
plt.legend(loc = 'best')


print('question 3 results graphed above')
   
    
