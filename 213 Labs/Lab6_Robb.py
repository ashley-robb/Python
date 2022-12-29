"""
ENPH 213 
Lab 6
Ashley Robb
"""
#%%

import matplotlib.pyplot as plt 
import numpy as np
import sympy as sm
from sympy import lambdify
from sympy.abc import h, t

#%%
'''
QUESTION 1

PART A
'''

def function_calculator(function ,t1 , t2, n): # creating a calculator to determine arrays of values
    if (n%2 != 0): # ending the function and ouputting an error if the n values is odd
        print('ERROR: the n value entered must be an even integer value (0,2,4,6,..etc.)') # printing error
        t = 'cannot provide values for the time output'
        y = 'cannot provide values for function output'
    else:
        t = np.arange(t1 , t2 , abs(t2-t1) / n) # creating t array 
        y = function(t) # plugging values from t array into function
    return t , y # returning time (x values) and y (function outputs)
     


# variables given for our use in this problem
a0, a1, a2 = 3, 1, 0.5
w0, w1, w2 = 1, 4, 7
    
#defining the given function
function1a = a0*sm.sin(w0*h) + a1*sm.sin(w1*h) + a2*sm.sin(w2*h)  # defining the function symbolically
f1a = lambdify([h], function1a) # lambdifying the function to allow us to plug in real numeric values

#determining t and y arrays to plot
q1a_30_t = function_calculator(f1a , 0 , 2*np.pi , 30)[0] # time values for n = 30
q1a_30_y = function_calculator(f1a , 0 , 2*np.pi, 30)[1] # corresponding function outputs for n = 30

q1a_60_t = function_calculator(f1a , 0 , 2*np.pi , 60)[0] # time vlaues for n = 60
q1a_60_y = function_calculator(f1a , 0 , 2*np.pi, 60)[1]  # corresponding function outputs for n = 60


def DFT(func, t1, t2, n): # defining DFT function 
    if (n%2 != 0): # stoping function if n is odd
        print('ERROR: the n value entered must be an even integer value (0,2,4,6,..etc.)') # error message
        y = 'cannot provide values for y array' 
        w = 'cannot provide omega values'
    else:
        x = function_calculator(func , t1, t2, n)[1] # getting x values from function calculator above
        k = np.linspace(0, n-1, num = n) # creating array of k values used for DFT evaluation
        i = np.linspace(0 , n-1, num = n) # creating an array of i values for DFT evaluation (called j in lecture slides)
        e = np.exp(-2j * np.pi * np.outer(k , i) / n) # defining the exponential 
        y = np.matmul(e, x) # matrix multiplication of values and exponents
        f = abs(i / (t2-t1)) # determining the frequency 
        w = (2 * np.pi * f) # getting omega using frequency 
    return y,w

def IDFT(x_vals, t1, t2, n): # creating IDFT function
    if (n%2 != 0): # stoping function if n is odd
        print('ERROR: the n value entered must be an even integer value (0,2,4,6,..etc.)') # error message
        y = 'cannot provide values for y array' 
        freq = 'cannot provide frequency values'
    else:
        x = x_vals # getting x values from function calculator above
        k = np.linspace(0, n-1, num = n) # creating array of k values used for IDFT evaluation
        i = np.linspace(0 , n-1, num = n) # creating an array of i values for IDFT evaluation (called j in lecture slides)
        e = np.exp(2j * np.pi * np.outer(k, i) / n) # defining the exponential  
        y = (np.matmul(e, x)) / n # matrix multiplication of values and exponents divided by n
        freq = np.arange(t1 , t2 , abs(t2-t1) / n) # determining the frequency 
    return y,freq

y_30 = DFT(f1a, 0 , 2*np.pi , 30)[0] # setting y values for n = 30 DFT
omega30 = DFT(f1a , 0 , 2*np.pi , 30)[1] # setting omega values for n=30 DFT

y_60 = DFT(f1a, 0 , 2*np.pi , 60)[0] # setting y vals for n = 60 DFT 
omega60 = DFT(f1a , 0 , 2*np.pi , 60)[1] # setting omega vals for n = 60 DFT


plt.subplot(1,2 ,1) # subplot format

plt.plot(q1a_60_t , q1a_60_y, 'b', label = 'n=60') # plot function for n = 60
plt.plot(q1a_30_t , q1a_30_y, 'r--', label = 'n=30') # plot function for n = 30 
plt.legend(loc = 'best') # adding legend
plt.ylabel('y', size =15) # axis labels and enlarged size
plt.xlabel('t', size = 15)
plt.text(-1, 5.2, 'QUESTION 1\nPART A:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') # question label 
plt.title('Graph of \nGiven Function')

plt.subplot(1, 2 ,2) #subplot format 

plt.stem(omega60, abs(y_60) , 'b' , markerfmt = 'none', basefmt = '-b' , label = 'n = 60') # stem function to show amplitudes only for n = 60 
plt.stem(omega30, abs(y_30) , 'r--', markerfmt = 'none', basefmt = '-r' ,  label = 'n = 30') # same for n = 30
plt.legend(loc = 'best') # add legend
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.title('DFT Grpah of \nGiven Function')
plt.show()

print('For part 1 a results, view the first plot')

#%%
'''
QUESTION 1

PART B
'''

'''
in part a I used the arange function to plot my graphs because it is more 
accurate than the linspace function. The linspace function adds small 
frequencies to each value in the array as it progresses because the function 
repeats itself. The arange function does not overlap the first and last value
the same way that the linspace function does. Included below is the different
plots provided based on whether linspace or arange was used to produce the 
plots. 
'''
# defining given vals 
t1 = 0
t2 = 2* np.pi 
n1 = 30
n2 = 60


# copied and pasted function from 1a changed time step generator 
def DFT_LS(func, t1, t2, n): # CHANGE: function name. defining DFT function 
    if (n%2 != 0): # stoping function if n is odd
        print('ERROR: the n value entered must be an even integer value (0,2,4,6,..etc.)') # error message
        y = 'cannot provide values for y array' 
        w = 'cannot provide omega values'
    else:
        x = np.linspace(t1, t2, n) #CHANGE: producing t values with linspace instead of calling from function calculator 
        k = np.linspace(0, n-1, num = n) # creating array of k values used for DFT evaluation
        i = np.linspace(0 , n-1, num = n) # creating an array of i values for DFT evaluation (called j in lecture slides)
        e = np.exp(-2j * np.pi * np.outer(k , i) / n) # defining the exponential 
        y = np.matmul(e, x) # matrix multiplication of values and exponents
        f = abs(i / (t2-t1)) # determining the frequency 
        w = (2 * np.pi * f) # getting omega using frequency 
    return y,w

y_30_LS = DFT_LS(f1a, 0 , 2*np.pi , 30)[0] # setting y values for n = 30 DFT with linspace
omega30_LS = DFT_LS(f1a , 0 , 2*np.pi , 30)[1] # setting omega values for n=30 DFT with linspace

y_60_LS = DFT_LS(f1a, 0 , 2*np.pi , 60)[0] # setting y vals for n = 60 DFT with linspace
omega60_LS = DFT_LS(f1a , 0 , 2*np.pi , 60)[1] # setting omega vals for n = 60 DFT with linspace


plt.subplot(1,2 ,1) # subplot format
# linspace vs arange for n = 30
plt.stem(omega30, abs(y_30) , 'b' , markerfmt = 'none', basefmt = '-b' , label = 'arange') # stem function to show amplitudes for n = 30 arange 
plt.stem(omega30_LS, abs(y_30_LS) , 'r--', markerfmt = 'none', basefmt = '-r' ,  label = 'linspace') # same for linspace
plt.title('   Linspace and Arange\n Function For n = 60')
plt.legend(loc = 'best') # add legend
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.text(-5, 105, 'QUESTION 1\nPART B:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') #question label 

plt.subplot(1, 2 ,2) #subplot format 
# linspace vs arange for n = 60
plt.stem(omega60, abs(y_60) , 'b' , markerfmt = 'none', basefmt = '-b' , label = 'arange') # stem function to show amplitudes for n = 60 arange
plt.stem(omega60_LS, abs(y_60_LS) , 'r--', markerfmt = 'none', basefmt = '-r' ,  label = 'linspace') # same for linspace
plt.title('Linspace and Arange\n Function For n = 60')
plt.legend(loc = 'best') # add legend
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('When viewing plots for part 1 b we can clearly see that the linspace function provides us with additional and incorrect frequencies')

#%%
'''
QUESTION 1

PART C
'''

idft_30 = IDFT(y_30 , 0, 2*np.pi, 30) # calling the inverse DFT function for n = 30
idft_60 = IDFT(y_60 , 0 , 2*np.pi, 60) # calling the inverse DFT function for n = 60 

plt.subplot(1,2 ,1) # subplot format

plt.plot(q1a_30_t , q1a_30_y, 'b', label = 'original function (n = 30)') # plot original function
plt.plot(idft_30[1] , idft_30[0], 'r-.', label = 'function transformed (n = 30)') # plot function through DFT and then back through IDFT
plt.legend(loc = 'best') # adding legend
plt.ylabel('y', size =15) # axis labels and enlarged size
plt.xlabel('t', size = 15)
plt.text(-2,5, 'QUESTION 1\nPART C:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') # question label 
plt.subplot(1, 2 ,2) #subplot format 

plt.plot(q1a_60_t , q1a_60_y, 'b', label = 'original function (n=60)') # plot original function
plt.plot(idft_60[1] , idft_60[0], 'r-.', label = 'function transformed (n=60)') # plot function through DFT and then back through IDFT
plt.legend(loc = 'best') # adding legend
plt.xlabel('t', size = 15)
plt.suptitle('Graphs of Given Function Comparing Original\nValues and Values Through DFT / IDFT Functions')
plt.show()

print('Notice for the plotted question 1 part c results, the original function and the function that has been transformed through the DFT ans back through the IDFT functions match up perfectly! yay :)')

#%%
'''
QUESTION 2

PART A
'''
#inputting given variables
sigma = 0.5
wp = 0
n = 60
t1 = -np.pi
t2 = np.pi

#defining the given function
function2 = sm.exp(-(t**2)/ (sigma**2)) * sm.cos(wp * t)  # defining the function symbolically
f2 = lambdify([t], function2) # lambdifying the function to allow us to plug in real numeric values


def gaus_pulse(t1 , t2, n):
    time = np.arange(t1 , t2, abs(t2-t1) / n)
    function2 = sm.exp(-(t**2)/ (sigma**2)) * sm.cos(wp * t)  # defining the function symbolically
    f2 = lambdify([t], function2) # lambdifying the function to allow us to plug in real numeric values
    return time , (f2(time))


pulse = gaus_pulse(-np.pi , np.pi , 60) # using new gaussian pulse function

plt.subplot(1,2,1) # subplot format 
plt.plot(pulse[0] , pulse[1], 'b')# plotting the function
plt.xlabel('t' , size = 15) # axis labels and enlarged font 
plt.ylabel('y(t)' , size = 15)
plt.title('Gaussian Function') # giving graph a title 
plt.text(-4,1.1, 'QUESTION 2\nPART A:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') # question label 

pulse_DFT = DFT(f2, -np.pi , np.pi, 60) # plugging values into DFT 

time = np.arange(t1 , t2, abs(t2-t1) / n) # defining times to solve dt
dt = time[1] - time[0] # using times to solve dt 

# modified notes code for shift 
w_shift = np.fft.fftfreq(n,dt)*2.*np.pi
w_shift = np.fft.fftshift(w_shift)
y_shift = np.fft.fftshift(pulse_DFT[0])


plt.subplot(1 , 2, 2) # subplot function 
plt.plot(pulse_DFT[1], abs(pulse_DFT[0]), 'b' , label = 'no shift') # plotting DFT original 
plt.plot(w_shift, abs(y_shift) , 'r--' , label = 'with shift') # plotting DFT shifted 
plt.legend(loc = 'best')
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.title('DFT of Gaussian Function\nWith and Without Shift')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('View Plots for results')

#%%
'''
QUESTION 2

PART B
'''
#inputting given variables
sigma = 1
wp1 = 10
wp2 = 20
n = 400

def gaus_pulse(t1 , t2, w, n): # defining a function for the gaussian pulse 
    time = np.arange(t1 , t2, abs(t2-t1) / n)
    function2 = sm.exp(-(t**2)/ (sigma**2)) * sm.cos(w * t)  # defining the function symbolically
    f2 = lambdify([t], function2) # lambdifying the function to allow us to plug in real numeric values
    return time , (f2(time)) # returning time array and function values 

pulse10 = gaus_pulse(-np.pi , np.pi , wp1, 400) # gaussian pulse for omega = 10 
pulse20 = gaus_pulse(-np.pi , np.pi , wp2, 400) # gaussian pulse for omega = 20 

plt.subplot(1,2,1) # subplot function 
plt.plot(pulse10[0] , pulse10[1], 'b', label = '$\omega_0$ = 10') # plotting the gaussian frequencies for omega = 10
plt.plot(pulse20[0] , pulse20[1] , 'r' , label = '$\omega_0$ =20') # same with omega = 20 
plt.xlabel('t' , size = 15) # labling axis and enlarging font 
plt.ylabel('y(t)' , size = 15)
plt.legend(loc = 'best')
plt.title('Gaussian Function at\n Different Frequencies') # title 
plt.text(-5,1.3, 'QUESTION 2\nPART B:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') # question box label 

function2 = sm.exp(-(t**2)/ (sigma**2)) * sm.cos(wp1 * t)  # redefining the function with omega = 10 
f2_1 = lambdify([t], function2) # lambdifying the redefined function 

pulse_DFT10 = DFT(f2_1, -np.pi , np.pi, n) # plugging values into DFT function 
time = np.arange(t1 , t2, abs(t2-t1) / n) # time array 
dt = time[1] - time[0] # solving dt 
#code modified from lectures 
w_shift = np.fft.fftfreq(n,dt)*2.*np.pi
w_shift10 = np.fft.fftshift(w_shift)
y_shift10 = np.fft.fftshift(pulse_DFT10[0])

function2 = sm.exp(-(t**2)/ (sigma**2)) * sm.cos(wp2 * t)  # redefining the function with omega = 20 
f2_2 = lambdify([t], function2) # lambdifying the redefined function
pulse_DFT20 = DFT(f2_2, -np.pi, np.pi, n) # plugging new omega = 20 vals into DFT
# code modified from lectures 
w_shift = np.fft.fftfreq(n,dt)*2.*np.pi
w_shift20 = np.fft.fftshift(w_shift)
y_shift20 = np.fft.fftshift(pulse_DFT20[0])

plt.subplot(1 , 2, 2) # subplot function 
plt.plot(pulse_DFT10[1], abs(pulse_DFT10[0]), 'b' , label = '$\omega_0$ = 10 (no shift)') # plotting omega = 10 without shift 
plt.plot(w_shift10, abs(y_shift10) , 'b' , linestyle='-', marker='.' , label = '$\omega_0$ = 10 (with shift)') # omega = 10 shifted 
plt.plot(pulse_DFT20[1], abs(pulse_DFT20[0]), 'r' , label = '$\omega_0$ = 20 (no shift)') # omega = 20 no shift
plt.plot(w_shift20, abs(y_shift20) , 'r' , linestyle='-', marker='.' , label = '$\omega_0$ = 20 (with shift)') # omega = 20 with shift 
plt.xlim([-40,40]) # setting x axis 
plt.legend(loc = 'best' , fontsize = 8)
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.title('DFT of Gaussian Function\nWith and Without Shift')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('View Plots for results')

#%%
'''
QUESTION 3
'''
# variables given for our use in this problem
a0, a1, a2 = 3, 1, 0
w0, w1, w2 = 1, 10,0

function3 = a0*sm.sin(w0*h) + a1*sm.sin(w1*h) + a2*sm.sin(w2*h)  # defining the function symbolically
f3 = lambdify([h], function3) # lambdifying the function to allow us to plug in real numeric values

  
q3_t = function_calculator(f3 , 0 , 8*np.pi , 200)[0] # time values for n = 200
q3_y = function_calculator(f3 , 0 , 8*np.pi, 200)[1] # corresponding function outputs for n = 200

y3 = DFT(f3, 0 , 8*np.pi, 200)[0] # question 3 DFT y vals 
omega3 = DFT(f3, 0 , 8*np.pi, 200)[1] # question 3 DFT omega vals 

y3_new = y3.copy() # copying the array 



for i in range(0, len(omega3)): # loop to eliminate noise 
    if omega3[i] > 5: # omega values after initial amplitude are noise 
        y3_new[i] = 0 # setting noise to zero 

y3_scaled = y3_new * 2 # y must be doubled to allow for amplitude of filtered IDFT to be correct
# the reason y is doubled is because the DFT scales everything relatively and there are 2 a0 values on our DFT graph 

idft_3 = IDFT(y3_scaled , 0, 8*np.pi, 200) # plugging values into IDFT 

plt.subplot(1,2 ,1) # subplot format
plt.plot(q3_t , q3_y, 'r--', label = 'unfiltered') # unfiltered function 
plt.plot(idft_3[1] , idft_3[0], 'b', label = 'filtered') # plot function through DFT, filtered, and then back through IDFT
plt.legend(loc = 'lower left') # adding legend
plt.ylabel('y', size =15) # axis labels and enlarged size
plt.xlabel('t', size = 15)
plt.text(-1, 5.2, 'QUESTION 3:', fontsize = 8, bbox=dict(facecolor='red') , horizontalalignment='center') # question label 
plt.title('Graph of \nGiven Function')

plt.subplot(1, 2 ,2) #subplot format 

plt.plot(omega3, abs(y3) , 'r--', label = 'unfiltered') # showing all DFT  frequencies  
plt.plot(omega3, abs(y3_new) , 'b',  label = 'filtered') # showing only DFT frequencies that arent noise 
plt.legend(loc = 'best') # add legend
plt.ylabel('|ỹ|' , size = 15) # axis labels with enlarged size
plt.xlabel('$\omega$' , size = 15)
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.title('DFT Grpah of \nGiven Function')
plt.show()

print('\n\n----------RESULTS FROM QUESTION 3: ----------\n')

print('The blue line in the graph provided is not scaled correctly. To scale it properly we need to multiply the DFT y values by 2 to accomoadate for the scaling the DFT function does and the frequencies we removed')

print('\n\nview outputted plots for all results')

print('\n\n----------All done! YAY :-) ----------')
