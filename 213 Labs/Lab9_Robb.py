#%% INFO
"""
ENPH 213 
LAB 9
NAME: ASHLEY ROBB 
STUDENT NUMBER: 20203465
DATE: 03/29/22
"""
#%%LIBRARIES
import matplotlib.pyplot as plt 
import numpy as np

#%% Q1
'''
QUESTION 1
'''
print('\n\n\n----------RESULTS FOR QUESTION 1: ----------\n\n')

n = 1000 # numper of points 
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, dpi = 300, figsize = (6,6)) # plotting layouts 
axis = np.array([ax1 , ax2, ax3, ax4]) # putting axis values into an array 
pi_vals = np.zeros(4) # array of four zeros to be replaced by pi values 


for i in range(0,4): # creating a loop to run through each of the subplots 
    rando = np.random.rand(n,2) # creating a 1000 x 2 array of random values to serve as x and y coordinates 
    rando_x = rando[:,0] - 0.5 # adjusting the values from (0,1) to (-.5, .5)
    rando_y = rando[:,1] - 0.5
    rando_radius = np.sqrt(rando_x**2 + rando_y**2) # determining the radial position of the point 
    circle = plt.Circle((0, 0), 0.5, color='b', alpha = 0.5) # creating the circle for the plot 
    count = 0 # starting a new count for each subplots pie value 
    for j in range(n): # running through all 1000 points 
        if rando_radius[j] < 0.5: # determining if the radial value of each of the points is within the circles radius 
            count = count + 1 # adding to the count if it is 
    pi_vals[i] = count * 4/n # storing the pi values for each subplot 
    plt.xlim(-0.6, 0.6) # plotting details 
    plt.ylim(-0.6,0.6)
    axis[i].title.set_text('pi = ${}$'.format(pi_vals[i])) # adding the pi value to the top of each of the plots
    axis[i].set_aspect(1) # formatting the axis 
    axis[i].add_patch(circle) # adding the cirle 
    axis[i].scatter(rando_x, rando_y, c=np.random.rand(len(rando_x),3), marker = '.') # scattering the points with random colours 
    print('pi value for plot',i+1, '=' , pi_vals[i]) # printing the pi values 
    
plt.text(-2,2.25, 'QUESTION 1' , fontsize = 12, bbox=dict(facecolor='red') , horizontalalignment='center') # question label
plt.show() # showing the subplots 



#%% Q2
'''
QUESTION 2
'''
print('\n\n\n----------RESULTS FOR QUESTION 2: ----------\n\n')

# defining the given values 
N = 20  
p = 0.6 
kbT1 = 0.1
KT1 = 1/kbT1

def initialize(N):
    spin = np.ones(N)
    E = 0; M = 0
    for i in range(1 , N ):
        if np.random.rand(1) < p :
            spin[i] = - 1
        E = E - spin[ i - 1 ] * spin[ i ] # Energy
        M = M + spin[i] # Magnetization
    # periodic bc inclusion
    E = E - spin [N - 1] * spin [0]
    M = M + spin[0]
    return spin , E , M

def update(N , spin , kT , E , M ):
    num = np.random.randint(0 , N - 1)
    flip = 0
    # periodic bc returns 0 if i + 1 == N , else no change :
    dE = 2 * spin[num]*(spin[num- 1] + spin[(num + 1) % N ])
    # if dE is negative , accept flip :
    if dE < 0 :
        flip = 1
    else :
        p = np.exp( - dE / kT )
        if np.random.rand(1) < p:
            flip = 1
    # otherwise , reject flip
    if flip == 1 :
        E = E + dE
        M = M - 2 * spin[num]
        spin[num] = -spin[num]
    return E,M

init = initialize(N)
spins = init[0]
E_vals = init[1]
M_vals = init[2]

new = update(N, spins, kbT1, E_vals, M_vals)

E1 = new[0]
M1 = new[1]

print(E1 , M1) # testing the function 

#%% 

print('NOTE: I attended the virtual labs in the first 6 weeks with my camera on! :)\nThanks!!')