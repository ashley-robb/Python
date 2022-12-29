"""
ENPH 213
ASHLEY ROBB
20203465
LAB 8
WINTER 2022
"""

#%% libraries 
import matplotlib.pyplot as plt 
import numpy as np
from numpy import linalg as lg
from timeit import default_timer

#%% QUESTION 1

'''
QUESTION 1
'''
print('\n\n\n----------RESULTS FOR QUESTION 1: ----------\n\n')

# declaring given values
t0 = 0# initial time
tf = 61 # final time 
dt = 0.1 # time step 

x0 = 0 # initial position
xf = 1 # final position
dx = 0.01 # position step 

t = np.arange(t0 , tf, dt) # array of time values
x = np.arange(x0, xf , dx) # array of position values 
alpha = 2.3 * 10**(-4) # setting alpha 
kappa = alpha * dt / (dx**2) # determining kappa with alpha 
plot_t = np.array([0,5,10,20,30,60]) # array of t values plotted



for n in t: # running through time array 
    if n == 0:
        u = 20 + 30 * np.exp(-100 * (x- 0.5)**2) #declaring initial u value
    else:
        u[1:-1] = u[1:-1] + kappa*(u[:-2] - 2*u[1:-1] + u[2:]) # other u values 
    if n in plot_t: # plotting the specific time values 
        plt.plot(x, u, label = '$time = {}$'.format(n))
                
plt.legend(loc = 'best') # adding legend 
plt.title('Heat Equation Solution at Various\n Snapshots in Time (Time Units are Seconds)', size = 14) # adding a title 
plt.xlabel('Rod length (m)', size = 14); plt.ylabel('Temperature (C)', size = 14) # axis labels 
plt.text(-0.1, 55, 'Figure 1', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
        
plt.show()
print('Refer to Figure 1 for question 1 results.')
#%% QUESTION 2

'''
QUESTION 2
'''
print('\n\n\n----------RESULTS FOR QUESTION 2: ----------\n\n')


x0 = 0 # initial x position 
xf = 2 # final x position 
dx = 1/50 # x position step 

y0 = 0 # initial y position 
yf = 1 # final y position 
dy = 1/50 # y position step 

x = np.arange(x0 , xf, dx) # array of x values
y = np.arange(y0, yf , dy) # array of y values

def func(x,y): # given function 
    return np.cos(10 * x) - np.sin(5 * y - np.pi/4)

def jacobi(xstep, ystep, xs, ys, function, tol): # jacobi function 
    tol_check = tol + 1 # setting the checker greater than the tolerance to gaurentee initial entry to while loop 
    x_pos,y_pos = np.meshgrid(xs,ys) # mesh grid from x and y position values 
    phi = np.zeros([50,100]) # creating a 2d array of zeros to be adjusted later
    fpq = function(x_pos, y_pos) # function at given points, needed for further evaluation 
    while(tol_check > tol): # checking if we are within tolerance 
        phi_save = np.copy(phi) # copying the grid 
        phi[1:-1, 1:-1] = (ystep**2*(phi[2:, 1:-1] + phi[:-2, 1:-1]) + xstep**2*(phi[1:-1, 2:] + phi[1:-1, :-2]) -(xstep**2 * ystep**2 *fpq[1:-1, 1:-1])) / (2 * (xstep**2 + ystep**2)) # given formula 
        tol_check = (lg.norm(phi) - lg.norm(phi_save))/ lg.norm(phi) # evaluating tolerance 
    return phi # returning output 

solution = jacobi(dx, dy, x, y, func, 10**(-5)) # calling function 


# plotting function with colour bar 
c = plt.imshow(solution, cmap ='hot', vmin = 0, vmax = 0.07,
                  extent =[x.min(), x.max(), y.min(), y.max()],
                    interpolation ='nearest', origin ='lower')
plt.colorbar(c)
                
plt.title("Solution to Poisson's Equation \nwith Fixed Boundary Conditions ", size = 14) # adding a title 
plt.xlabel('x', size = 14); plt.ylabel('y', size = 14) # axis labels 
plt.text(-0.1 , 1.2, 'Figure 2', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
     

plt.show()

print('Refer to Figure 2 for question 2 results.')
#%% QUESTION 3

'''
QUESTION 3
'''
print('\n\n\n----------RESULTS FOR QUESTION 3: ----------\n\n')

# declaring given values 
step = np.pi/ 400 # step size
x = np.arange(0 , 2* np.pi, np.pi/400) # x values     
y = np.arange(0 , 2* np.pi, np.pi/400)  # y values 
n = 800 # number of values
l = np.arange(n) # array 
k = np.arange(n) # array 

x_grid, y_grid = np.meshgrid(x,y) # declaring meshgrids 
k_grid , l_grid = np.meshgrid(k , l)


start_time = default_timer()#starting timer
# taking the fft of f(x,y)
def func_3(x,y): # given function 
    return np.cos(3 * x + 4 * y) - np.cos(5 * x - 2 * y)

func3 = func_3(x_grid , y_grid)
fft= np.fft.fft2(func3)

# plugging into phi tilde 
def phi_tilde(h, k, l, n, f_tilde):# defining the given function 
    phi_ti = np.zeros([n,n])
    phi_ti[1:,1:] = np.real((h**2 *f_tilde[1:,1:]) / (2 * (np.cos((2*np.pi *k[1:, 1:]) / n) + np.cos((2* np.pi * l[1:,1:]) / n) - 2))) # denominator 
    return phi_ti

phi_t = phi_tilde(step, k_grid, l_grid , n, fft)

# plugging phi tilde into ifft
ifft = np.fft.ifft2(phi_t)

stop_time = default_timer() # stopping timer because calculations are done 

print('The time this calculation takes is:' , round(stop_time - start_time, 4) , 'seconds\n') # printing results 

# plotting results 
q3 = plt.imshow(np.real(ifft), cmap ='hot', vmin = -.07, vmax = 0.07,
                  extent =[0,2,0,2],
                    interpolation ='nearest', origin ='lower')
bar = plt.colorbar(q3)
bar.set_label('$\phi(x,y)$' , size = 12)
plt.title("Solution to Poisson's Equation \nOn a Periodic Grid and Solved in a Fourier Space", size = 14) # adding a title 
plt.xlabel('$x(\pi)$', size = 14); plt.ylabel('y$(\pi)$', size = 14) # axis labels 
plt.text(-0.3 , 2.3, 'Figure 3', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label

plt.show()

print('Refer to Figure 3 for question 3 results.\n')

print('All done :) Thanks!')
























