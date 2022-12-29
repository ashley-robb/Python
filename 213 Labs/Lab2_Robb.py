#%%
'''
ASHLEY ROBB 
STUDENT ID: 20203465
213 LAB 2
DUE: 01/30/2022
'''
#%%
'''
QUESTION 1
PART A
'''
import numpy as np 
import math as m
from scipy import integrate 

print("---------------RESULTS FROM QUESTION 1 PART A:---------------\n")
#code given
n=100
a=0
b=1
h=(b-a)/(n-1)
int = 0 # Initialize the integral variable
for i in range(0,n-1):# fixing the range from (0,n) to (0,n-1)
    xi = a + h*i # Determine the xi value for the loop
    int = int + (2/np.pi **0.5)*np.exp(-xi **2)*h

print ('inefficient rectangle method:', int) # printing the 'bad code' results 
'''
NOTE: the following code was done for QUESTION 1 PART B, however it needed to 
be defined before the rect, trap, and simps functions so it has been moved here
'''

 
def reler_any(g,f,a,b,n):
    return abs((g-(integrate.quad(f,a,b)[0]))/(integrate.quad(f,a,b)[0]))


'''
RESUME QUESION 1 PART A 
'''   
def error(x): return (2/np.pi **0.5)*np.exp(-x **2) # defining our function for easy access

def rect(f,a,b,n): #defining a more efficient rectangle integration function 
    h=(b-a)/(n-1)# definfng h 
    ans = 0 # Initialize the integral variable
    for i in range(0,n-1): # setting the range for the loop
        x = a + h*i # Determine the x value for the loop
        ans = ans + f(x) *h # adding to the answer by calling the function 
    if (reler_any(ans, f, a, b, n) > 0.001): # creating an if statement to display an error if code isnt accurate 
        output = ans, 'WARNING: results have a relative error greater than 0.1%'
    else: output = ans, 'results accurate withing 0.1% relative error'
    return output

print('more efficient rectangle method:', rect(error, 0, 1, 100)) # printing results from effiecient rectangle integration function 

#%%
'''
QUESTION 1
PART B
'''  
print('\n\n---------------RESULTS FROM QUESTION 1 PART B:---------------\n')

def trap(f,a,b,n): #defining an efficient trapezoidal integration function 
    ans = 0. # Initialize the integral variable
    h=(b-a)/(n-1)
    for i in range(0,n-1): # setting the range for the loop
        xi = a + h*i # Determine the xi and xj value for the loop
        xj = a + h*(i + 1) 
        ans = ans + (f(xi) + f(xj)) *h*(1/2) #adding to the answer by calling the function 
    if (reler_any(ans, f, a, b, n) > 0.001): # using an if statement to display a warning for big errors
        output = ans, 'WARNING: results have a relative error greater than 0.1%'
    else: output = ans, 'results accurate within 0.1% relative error'
    return output


print('trap(n=100):', trap(error, 0,1, 100)) # printing trapezoidal method for n = 100
print('trap (n=101):', trap(error, 0, 1, 101)) # and for n = 101


def simps(f,a,b,n): #defining an efficient simpsons method integration function
    if (n%2 == 0): # giving an error for even n values 
        output = 'ERROR: simpsons rule needs an even number of panels and therefore n must be odd'
    else:
        ans = 0 # Initialize the integral variable
        h = (b-a)/(n-1) # definfing h
        for i in range(0, n-2, 2): # setting the range for the loop
            first = a + h*i # first value in equation 
            mid = a +h*(i+1) #second value in equation
            last = a + h*(i+2)# last value in equation 
            ans = ans + h*(1/3)*(f(first)+4*f(mid)+f(last)) # adding to the answer by calling the function 
        if (reler_any(ans, f, a, b, n) > 0.001): # printing a warning for relative error larger than .1%
            output = ans, 'WARNING: results have a relative error greater than 0.1%'
        else: output = ans, 'results accurate withing 0.1% relative error'
    return output

print('simpsons method (n=100):', simps(error, 0, 1, 100)) # printing simpsons method for n = 100 (notice the error message)
print('simpsons method (n=101):', simps(error, 0, 1, 101)) # printing simpsons for n = 101

print('python math error function:', m.erf(1)) # printing the python function to compare

def reler_erf(g,f,a,b,n): #defining a relative error function for our function with erf
    return abs((g(f,a,b,n)[0]-m.erf(b-a))/m.erf(b-a))*100


print('rectangle method relative error for n = 100:', (reler_erf(rect,error,0,1,100)), '%') # printing relative error 
print('rectangle method relative error for n = 101:', (reler_erf(rect,error,0,1,101)), '%')

print('trapezoidal method relative error for n = 100:', (reler_erf(trap,error,0,1,100)), '%')
print('trapzoidal method realtive error for n = 101:', (reler_erf(trap,error,0,1,101)), '%')

#simpson function does not work for n = 100
print('simpsons method realtive error for n = 101:', (reler_erf(simps,error,0,1,101)), '%')

#%%
'''
QUESTION 1
PART C
'''  
print('\n\n---------------RESULTS FROM QUESTION 1 PART C:---------------\n')
print('be patient question 1 c) takes about 10 seconds to output')
print('(see comment in line 145)')

def adaptive_step(k,f,a,b): # defining an adaptive step function to check relative error 
    ni = 3 # initial n 
    nj = 2*ni - 1 # n prime 
    av = abs(k(f,a,b,ni)[0]-k(f,a,b,nj)[0]) # absolute difference between function at n and n prime 
    if (k == trap): 
        E = (1/3) # setting the error fraction for trapezoidal 
    elif (k == simps):
        E = (1/15) # error fraction for simpsons 
    else:
        nj = 'ERROR: this function is only designed for trapezoidal and simpsons method'
        #returning an error for functions other than trap and simps
    while ( E*(av) > 1e-13): # while loop for unacceptable relative error 
        ni = 2*ni-1 # updating n 
        nj = 2*ni - 1 # updating n prime 
        av = abs(k(f,a,b,ni)[0]-k(f,a,b,nj)[0]) # updating av for the while loop 
    return nj # returning the number of n values required to have a function accurate to 10^-13

print('minimum numerical cell value for accuracy within 10^-13 relative error for trapezoidal function:', adaptive_step(trap, error, 0, 1))

print('minimum numerical cell value for accuracy within 10^-13 relative error for simpsons function:', adaptive_step(simps, error, 0, 1))

#The adaptive_step(trap, error, 0, 1) takes about 10 seconds to run while the rest of the code 
#takes a fraction of a second. If you know why this is please leave some feedback explaining :) Thanks !
#%%
'''
QUESTION 2
PART A
'''  
print('\n\n---------------RESULTS FROM QUESTION 2 PART A:---------------\n')
import pandas as pd # used to read file 
import matplotlib.pyplot as plt # used to plot data from file 

df = pd.read_csv(r"C:\Users\Ashle\Downloads\Hysteresis-Data.csv") # reads the data in the file 

vx = np.array(df['vx']) # converting vx and vy from pandas.core.series.Series to numpy.ndarray
vy = np.array(df['vy']) # the conversion was necessart for question 2 part b

plt.plot(vx,vy) # plotting the data 
plt.xlabel('$V_x$') # x axis label 
plt.ylabel('$V_y$') #y axis label 
plt.title('Hysteresis Data') # plot title 

print('view plot titled "Hysteresis Data"')

#%%
'''
QUESTION 2
PART B
'''  
print('\n\n---------------RESULTS FROM QUESTION 2 PART B:---------------\n')

def int_cl(a,b): # defining a function for integral of closed line 
    area = -(a[1:]-a[0:-1])*(b[1:]+b[0:-1])/2 # using trapezoidal method to determine area of slices 
    return np.sum(area) # summing the slices 
print('area inside the Hysteresis curve:' , int_cl(vx,vy)) # printing with given data 

#%%
'''
QUESTION 3
PART A
'''
print('\n\n---------------RESULTS FROM QUESTION 3 PART A:---------------\n')

def simps_2d(f_2d , a , b , c , d , n , m): #defining a function for 2d functions using simpsons method 
    if (n%2 == 0 or m%2 == 0): # giving an error if n or m are inputted as even numbers
        output = 'ERROR: simpsons rule needs an even number of panels and therefore n and m must be odd'
    else: # if n and m are odd code will run 
            h = (b-a)/(n-1) # h for x
            g = (d-c)/(m-1) # h for y 
            x_vals = np.linspace(a, b, n) # assigning n number of values evenly between a and b 
            y_vals = np.linspace(c, d, m) # same as above but with m, c and d
            x,y = np.meshgrid(x_vals, y_vals, indexing = 'ij') # turning the 1d arrays into a 2d grid 
            x_vector = np.zeros(n) # creating an array to store the weight matrix values for x
            y_vector = np.zeros(m) # same for y 
            x_vector[0::2] = h*(2/3); y_vector[0::2] = g*(2/3) # setting every other place in the array to (h/3)*2 
            x_vector[1::2] = h*(4/3); y_vector[1::2] = g*(4/3) # same with (h/3)*2
            x_vector[0] = h*(1/3); y_vector[0] = g*(1/3) # setting first values to h/3
            x_vector[-1] = h*(1/3); y_vector[-1] = g*(1/3) # setting last values to h/3
            matrix = np.outer(x_vector, y_vector) # computing the weight matrix 
            output = np.sum((f_2d(x,y)*matrix)) # summing values and multiplying by weighted matrix  
    return output

print('no output for question 3 part A')

#%%
'''
QUESTION 3
PART B
'''
print('\n\n---------------RESULTS FROM QUESTION 3 PART B:---------------\n')

def given_2df(x,y): # defining the function given for question 3 
    return (x**2 + y)**(1/2) * np.sin(x) * np.cos(y)
# printing simpsons results for the 2d integral with varying values of n and m
print('Simpsons 2D results (n, m =101, 101)', simps_2d(given_2df, 0, np.pi, 0, np.pi/2, 101, 101))

print('Simpsons 2D results (n, m =1001, 1001)', simps_2d(given_2df, 0, np.pi, 0, np.pi/2, 1001, 1001))

print('Simpsons 2D results (n, m = 51, 101)', simps_2d(given_2df, 0, np.pi, 0, np.pi/2, 51, 101))

def compare(A, B): # defining a function to compare my value to the value given 
    if A == B:
        output = 'results from my code match given results from Dr. S. Hughes'
    else:
        diff = abs(A-B)/A
        output = 'results from my code do NOT match Dr. S. Hughes results.The reletive error is:' , diff
    return output

print('comparing 2d simpsons for n, m = 1001, 1001', compare(3.5389940350753895, (simps_2d(given_2df, 0, np.pi, 0, np.pi/2, 1001, 1001))))
    
print('comparing 2d simpsons for n, m = 51,101',compare(3.5389937460658536, (simps_2d(given_2df, 0, np.pi, 0, np.pi/2, 51, 101))))
'''
 my results for simpsons 2d method match for 51,101 but for 1001,1001 my value differd by a very small amount
 I believe that the difference comes from the storing and rounding done by the computer 
 I can see that for larger values of n and m the results become more accurate, however the accuracy only 
improves noticibly if n and m are both large value 
(if only one is large, the function is only about as accurate as the smallest value)
'''
      
#%%
'''
QUESTION 3
PART C
'''
from mpmath import *

print('\n\n---------------RESULTS FROM QUESTION 3 PART C:---------------\n')
# defining a lambda function for the given equation in Q3
func_Q3 = lambda x,y: sqrt(x**2 + y ) * sin(x) * cos(y) 

print('value computed using lambda quad function' , quad(func_Q3, [0, mp.pi], [0, mp.pi/2])) # checking the answers with quad integration 

# my value is the same as Dr. S. Hughes 

#%%
'''
QUESTION 3
PART D
'''
print('\n\n---------------RESULTS FROM QUESTION 3 PART D:---------------\n')

print('value and error computed using dblquad function', integrate.dblquad(func_Q3, 0, np.pi/2, 0, np.pi)) # checking answers using the dblquad funciton 
# my value and error is the same as Dr. S. Hughes 




















