"""
ENPH 213 
LAB 7 
ASHLEY ROBB 
"""

#%%
import matplotlib.pyplot as plt 
import numpy as np
import math as m 
import sympy as sm
import scipy as sc
from scipy.optimize import fsolve 
import matplotlib.colors as mcolors

#%% Q1 A

'''
QUESTION 1
PART A
'''
print('\n\n\n----------RESULTS FOR QUESTION 1 PART A: ----------\n\n')


plt.rcParams.update({'font.size': 18}) # keep those graph fonts readable!
plt.rcParams['figure.dpi'] = 120 # plot resolution

# similar to Matlab's fval function - allows one to pass a function
def feval(funcName, *args):
    return eval(funcName)(*args)

# vectorized forward Euler with 1d numpy arrays
def eulerf(f, y0, t, h): # Vectorized forward Euler (so no need to loop) 
    k1 = h*f(y0,t)                     
    y1=y0+k1
    return y1 

# vectorized back Euler with 1d numpy arrays
def eulerb(f, y0, t, h): # Vectorized backward Euler (so no need to loop) 
    eulerb_y = lambda y: y - y0 - h*f(y , t+h)
    ys = fsolve(eulerb_y , t)
    return ys

# stepper function for integrating ODE solver over some time array
def odestepper(odesolver, deriv, y0, t):
# simple np array  
    y0 = np.asarray(y0) # convret just in case a numpy was not sent
    y = np.zeros((t.size, y0.size))
    y[0,:] = y0; h=t[1]-t[0]
    y_next = y0 # initial conditions 

    for i in range(1, len(t)):
        y_next = feval(odesolver, deriv, y_next, t[i-1], h)
        y[i,:] = y_next
    return y

def deriv(y,t):
    return -10*y # derivative of function

def fun(t):
    return np.exp(-10*t) # given function 

def plotme(y1,y2,y3,n,ts): # plot function given 
    plt.plot(ts, y1, '-r', label='Exact', linewidth=3)
    plt.plot(ts, y2, 'gs', label='F-Euler $n={}$'.format(n), markersize=6)
    plt.plot(ts, y3, 'bo', label='B-Euler $n={}$'.format(n), markersize=6)
    plt.xlabel('$t$')     
    plt.legend(loc='best')

#defining given values 
n1=10 
ts1=np.linspace(0,0.6,n1)
y0= 1
yf1 = odestepper('eulerf',deriv, y0, ts1)
yb1 = odestepper('eulerb',deriv, y0, ts1)
yexact1 = fun(ts1)

# plotting results 
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plotme(yexact1 ,yf1 ,yb1, n1, ts1)
plt.ylabel('$y$')
plt.text(-0.15, 1.3, 'Figure 1', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.text(0.03, 1.1, '(a)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label

# defining the given variables 
n2=20
ts2=np.linspace(0,0.6,n2)
yf2 = odestepper('eulerf',deriv, y0, ts2)
yb2 = odestepper('eulerb',deriv, y0, ts2)
yexact2 = fun(ts2)

# plotting the results 
plt.subplot(1,2,2)
plotme(yexact2, yf2, yb2 ,n2 , ts2)
plt.text(0.03, 1.1, '(b)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label

plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('For Question 1 A Results Refer to Figure 1:\n')
# results explanation
print('Figure 1: Forward and Backward Euler ODE solvers plotted with the exact function y(t) = e^(-10t) spanning a time array of t = [0, 0.6]')
print('    Figure 1 (a): n = 10 uniformly spaced points.')
print('    Figure 1 (b): n = 20 uniformly spaced points.')


#%% Q1 B

'''
QUESTION 1
PART B
'''


print('\n\n\n----------RESULTS FOR QUESTION 1 PART B: ----------\n\n')

def deriv2(t , initial_vals): # system of equations 
    dy = np.zeros(len(initial_vals))
    dy[0] = initial_vals[1]
    dy[1] = -initial_vals[0]
    return dy

def fun2(t):
    return np.cos(t) # given function 
# vectorized RK4 with 1d numpy arrays
def RK4(f, y0, t, h): # Vectorized RK4 function  (so no need to loop)   
    k0 = h*f(t,y0)         
    k1 = h*f(t+h/2,y0+k0/2)
    k2 = h*f(t+h/2,y0+k1/2)      
    k3 = h*f(t+h,y0+k2)
    y = y0 + 1/6*(k0+2*k1+2*k2+k3)
    return y

def eulerf_new(f, y0, t, h): # Vectorized forward Euler (so no need to loop) 
    k1 = h*f(t,y0)                     
    y1=y0+k1
    return y1 

def t(T0): # converting time to period 
    period = 2 * np.pi 
    return T0/(period)

initial_vals = np.array([1 , 0]) # setting given values 
dt1 = 0.01
dt2 = 0.005
t1 = np.arange(0, 10* (2*np.pi), dt1) # time arrays 
t2 = np.arange(0, 10 *(2*np.pi), dt2)

#plotting results FOR DT = 0.01
rk4_1 = odestepper('RK4', deriv2, initial_vals, t1)
fe_1 = odestepper('eulerf_new',deriv2, initial_vals, t1)

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(rk4_1[:,0] , rk4_1[:,1], 'b', label = 'RK4')
plt.plot(fe_1[:,0] , fe_1[:,1] , 'r--', label = 'F-Euler')
plt.legend(loc = 'best')
plt.text(-2.3, 1.7, 'Figure 2', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.text(-1.13, 1.27, '(a)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('v'); plt.ylabel('x'); plt.title('dt = 0.01')


plt.subplot(1,2,2)
plt.plot(t(t1), rk4_1[:,0] ,'b', label = 'RK4' )
plt.plot(t(t1), fe_1[:,0], 'r--', label = 'F-Euler')
plt.plot(t(t1) , fun2(t1), '--' , color = 'yellowgreen', label = 'Exact')
plt.legend(loc = 'best')
plt.text(0.8, 1.28, '(b)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('t($T_0$)'); plt.ylabel('x')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('For Question 1 B Results Refer to Figure 2 and Figure 3:\n')
# explaning outputs 
print('Figure 2: 10 periods of simple harmonic oscillations solved with Fourth-Order Runge-Kutta and Forwards Euler methods, using a time step of dt = 0.01.')
print('    Figure 2 (a): Phase space solution to the simple harmnic ocillator.')
print('    Figure 2 (b): Position vs. Time plot of the same simple harmonic oscillator.')
#plotting results FOR DT = 0.005

rk4_2 = odestepper('RK4', deriv2, initial_vals, t2)
fe_2 = odestepper('eulerf_new',deriv2, initial_vals, t2)

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(rk4_2[:,1] , rk4_2[:,0], 'b', label = 'RK4')
plt.plot(fe_2[:,1] , fe_2[:,0] , 'r--', label = 'F-Euler')
plt.legend(loc = 'best')
plt.text(-2.3, 1.7, 'Figure 3', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.text(-.96, 1.1, '(a)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('v'); plt.ylabel('x'); plt.title('dt = 0.005')


plt.subplot(1,2,2)
plt.plot(t(t2), rk4_2[:,0] ,'b', label = 'RK4' )
plt.plot(t(t2), fe_2[:,0], 'r--', label = 'F-Euler')
plt.plot(t(t2) , fun2(t2), '--' , color = 'yellowgreen', label = 'Exact')
plt.legend(loc = 'best')
plt.text(0.8, 1.2, '(b)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('t($T_0$)'); plt.ylabel('x')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

# explaining outputs 
print('\nFigure 3: 10 periods of simple harmonic oscillations solved with Fourth-Order Runge-Kutta and Forwards Euler methods, using a time step of dt = 0.005.')
print('    Figure 3 (a): Phase space solution to the simple harmnic ocillator.')
print('    Figure 3 (b): Position vs. Time plot of the same simple harmonic oscillator.')
#%% Q2 A

'''
QUESTION 2
PART A
'''
print('\n\n\n----------RESULTS FOR QUESTION 2 PART A: ----------\n\n')
# defining given variables 
alpha_a = 0.0
beta_a = 1
gamma_a = 0.04
omega_a = 1
t0_a = 0; tf_a = 40*2*np.pi; dt_a = 0.01
F_a = 0.2
# setting array for y variabls 
y_a = np.zeros((2) , float)
y_a[0] = -.1; y_a[1] = 0.1

#y = np.array([-0.1, 0.1]) # initial x, initial v 
t_duff_a = np.arange(t0_a, tf_a, dt_a)

def func_a(t, y): # defining the given funciton 
    sln = np.zeros(len(y))
    sln[0] = y[1]
    sln[1] = -2 * gamma_a *y[1] - alpha_a *y[0] - beta_a*y[0]**3 + F_a*np.cos(omega_a*t)
    return sln

duff_a = odestepper('RK4', func_a, y_a, t_duff_a) # saving results 

duff_pos_a = duff_a[:,0] # position of results 
duff_vel_a = duff_a[:,1] # velocity of results 

# plotting results 
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(duff_vel_a[len(duff_vel_a)//4:] , duff_pos_a[len(duff_pos_a)//4:], 'r')
plt.text(-0.5, 0.55, 'Figure 4', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.text(-.33, 0.45, '(a)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.plot(duff_vel_a[0] , duff_pos_a[0], 'bo', markersize = 12)
plt.plot(duff_vel_a[-1::] , duff_pos_a[-1::], 'go' , markersize = 12)
plt.title('D Osc: $\\alpha$, F, $\omega$ = {},{},{}'.format(alpha_a, F_a, omega_a))
plt.xlabel('v'); plt.ylabel('x')


plt.subplot(1,2,2)
plt.plot(t(t_duff_a), duff_a[:,0] ,'b', linewidth = 1)
plt.text(0.2, 0.49, '(b)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('t($T_0$)'); plt.ylabel('x')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('F = 0.2, for results see Figure 4 plot')

# explaining output plots 
print('\nFigure 4: Duffing Oscillator')
print('    Figure 4 (a): Phase space solution to the duffing ocillator.')
print('    Figure 4 (b): Position vs. Time plot of the same duffing oscillator.')

#%% Q2 B

'''
QUESTION 2
PART B
'''
print('\n\n\n----------RESULTS FOR QUESTION 2 PART B: ----------\n\n')

# setting the given varaibles (others are reused from 2a, no need to redeclare)
alpha_b = 0.1
F_b = 7.5

def func_b(t, y): # defining the given function 
    sln = np.zeros(len(y))
    sln[0] = y[1]
    sln[1] = -2 * gamma_a *y[1] - alpha_b *y[0] - beta_a*y[0]**3 + F_b*np.cos(omega_a*t)
    return sln

duff_b = odestepper('RK4', func_b, y_a, t_duff_a) # save results 

duff_pos_b = duff_b[:,0] #results position 
duff_vel_b = duff_b[:,1] # results velocity 

# plotting results 
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(duff_vel_b[len(duff_vel_b)//2:], duff_pos_b[len(duff_pos_b)//2:], 'r')
plt.text(-8.5, 5, 'Figure 5', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.text(-6.3, 3.4, '(a)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.plot(duff_vel_b[0] , duff_pos_b[0], 'bo', markersize = 12)
plt.plot(duff_vel_b[-1::] , duff_pos_b[-1::], 'go' , markersize = 12)
plt.title('D Osc: $\\alpha$, F, $\omega$ = {},{},{}'.format(alpha_b, F_b, omega_a))
plt.xlabel('v'); plt.ylabel('x')


plt.subplot(1,2,2)
plt.plot(t(t_duff_a), duff_b[:,0] ,'b', linewidth = 1)
plt.text(0.2, 3.4, '(b)', fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.xlabel('t($T_0$)'); plt.ylabel('x')
plt.tight_layout() # tight layout is used for subplotted graphs to allow them to be as close as possible without the y axis label of the second plot overlapping with the first plot
plt.show()

print('F = 7.5, for results see Figure 5 plot')
# explaining output plots 
print('\nFigure 5: Duffing Oscillator')
print('    Figure 5 (a): Phase space solution to the chaotic duffing ocillator.')
print('    Figure 5 (b): Position vs. Time plot of the same chaotic duffing oscillator.')

#%% Q3 A

'''
QUESTION 3
PART A
'''
print('\n\n\n----------RESULTS FOR QUESTION 3 PART A: ----------\n\n')
# defining the known variables 
gamma = 10
r = 28
b = 8/3
t0 = 0
tf = 2*4*np.pi
dt = 0.01

# setting vectorized values 
time = np.arange(t0,tf,dt)
vals_1 = np.array([10,10,10])
vals_2 = np.array([1,1,1])

def func3a(t, vals): # system of equations 
    x = np.zeros(len(vals))
    x[0] = gamma * (vals[1] - vals[0])
    x[1] =  r *vals[0] - vals[1] - vals[2]*vals[0]
    x[2] = vals[0] * vals[1] - b*vals[2]
    return x
# saving values 
y1 = odestepper('RK4', func3a, vals_1, time)
y2 = odestepper('RK4', func3a, vals_2, time)


# plotting outputs 
ax = plt.figure(dpi = 200).add_subplot(projection = '3d')
ax.view_init(azim = 20 , elev = 29)
ax.plot(y1[:,0], y1[:,1], y1[:,2], linewidth = 0.5, label = 'y0 = [10,10,10]')
ax.plot(y2[:,0], y2[:,1], y2[:,2], linewidth = 0.5, label = 'y0 = [1,1,1]')
ax.set_title('Lorentz Attractor')
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
ax.legend(loc = 'best', fontsize = 8)
plt.show()

#%% Q3 B

'''
QUESTION 3
PART B
'''
print('\n\n\n----------RESULTS FOR QUESTION 3 PART B: ----------\n\n')
# importing new libraries 
from matplotlib import animation
#from IPython.display import HTML
import numpy as np
import matplotlib.pyplot as plt

#REST OF CODE GIVEN/ SLIGHTLY MODIFIED
# rename x,y,z for clarity in settign up simulations
x11=np.array(y1[:,0]) 
y11=np.array(y1[:,1])
z11=np.array(y1[:,2])
x22=np.array(y2[:,0]) 
y22=np.array(y2[:,1])
z22=np.array(y2[:,2])


time_vals = time
plt.rcParams.update({'font.size': 18})
fig = plt.figure(dpi=180) # nice and sharp!
ax = fig.add_axes([0.1, 0.1, 0.85, 0.85], projection='3d')
line, = ax.plot3D(x11, y11, z11, 'r-', linewidth=0.8)
line2, = ax.plot3D(x22, y22, z22, 'b-', linewidth=0.8)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

def init():
    line.set_data(np.array([]), np.array([]))
    line.set_3d_properties([])
    line.axes.axis([-2, 2, -2, 2])
    line2.set_data(np.array([]), np.array([]))
    line2.set_3d_properties([])
    line2.axes.axis([-2, 2, -2, 2])   
    return line,line2

def update(num):
    line.set_data(x11[:num], y11[:num])
    line.set_3d_properties(z11[:num])
    line2.set_data(x22[:num], y22[:num])
    line2.set_3d_properties(z22[:num])
    fig.canvas.draw()
    return line,line2

ani = animation.FuncAnimation(fig, update, init_func=init, interval=1, frames = len(time_vals), blit=True, repeat=True)


"""Uncomment if you want a movie, but need to wait a few minutes until finished, before the
screen one starts"""
# # conda install -c conda-forge ffmpeg # if you want to make an mp4
# # repeat Ture is teh default
# Writer = animation.writers['ffmpeg']
# #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=180)
# writer = Writer(fps=30)
# ani.save('im.mp4', writer=writer)

plt.show()







































