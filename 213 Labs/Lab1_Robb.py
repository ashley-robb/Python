#%%
'''
QUESTION 1
PART A
'''


#%% reminder to myself for how a matrix is numbered 

print('(0,0)    (0,1)    (0,2)    (0,3) \n(1,0)    (1,1)    (1,2)    (1,3)\n(2,0)    (2,1)    (2,2)    (2,3)\n(3,0)    (3,1)    (3,2)    (3,3)\n')


#%%
import sympy as sm # using sympy instructed 
from sympy import * 
from sympy.abc import x # using x as a variable 
import cmath as m
import numpy as np
import matplotlib.pyplot as plt
import math as ma

label = ['given function: ', 'first derivative: ', 'second derivative: ', \
         'third derivative: '] # labels each derivative when printing 
# learned to split long lines from:  https://itsmycode.com/syntaxerror-unexpected-character-after-line-continuation-character/

derivative = lambda i: (sm.diff(sm.exp(sm.sin(2*x)), x, i)) # creating a lamda function for the derivtive of the given function 
                  
for i in range(0,4): # creating a loop for calculating and printing derivative 
    print(label[i], derivative(i)) # prints the labels and their respective derivatives
    
f = lambdify(x,derivative(0),'numpy') # 23-26 converting each derivative from sympy to numpy
fd = lambdify(x,derivative(1),'numpy') #fd - first derivative
sd = lambdify(x,derivative(2),'numpy') # sd- second derivative
td = lambdify(x,derivative(3),'numpy') # td- third derivative 

a = np.linspace(0, 2*m.pi, 200) #creating an array with 200 elements evenly spaced between 0 and 2 pi


print(len(a), a[0], a[len(a)-1]-2*m.pi) # print test given in assignment 

function = np.array(f(a)) #33-36 creating arrays for the function and each derivaitve and evaluating them from 0 to 2 pi
first = np.array(fd(a))
second = np.array(sd(a))
third = np.array(td(a))

plt.plot(a,function, label='f(x)') #38-41 plotting the function and derivatives, as well as declaring labels  
plt.plot(a,first, label = "f'(x)")
plt.plot(a,second, label = "f''(x)")
plt.plot(a,third, label = "f'''(x)")

plt.xlabel('function input') # x axis label 
plt.ylabel('exp[sin(2x)] and its first three derivatives') #y axis label 
plt.title('f(x) = exp[sin(2x)] and its first three derivatives plotted from 0 to 2pi') # graph title
plt.legend() # adding a legend with the labels declared in 38-41
plt.show()# showing the graph 

#%%
'''
QUESTION 1
PART B 
'''
def centd(x,h): # defining a function for central differencing 
    centdiff = ((f(x+(h/2)) - f(x-(h/2))) / h) 
    return centdiff

def ford(x,h): #defining a function for forwards differencing 
    fordiff = ((f(x+h) - f(x)) / h) 
    return fordiff

a = np.linspace(0, 2*m.pi, 200) #creating an array with 200 elements evenly spaced between 0 and 2 pi

fordifs = np.array(ford(a, 0.15)) # first order forward differencing derivative- small step
centdifs = np.array(centd(a, .15))  # first order central differencing derivative - small step  

fordifb = np.array(ford(a,.5)) # first order forward differencing derivative - big step
centdifb = np.array(centd(a, .5))  # first order central differencing derivative - big step 


plt.plot(a,fordifs, label='forward differencing (.15)') #71-73 plotting the forward differencing and the analytical solution  
plt.plot(a,fordifb, label = "forward differencing (.5)") 
plt.plot(a,first,'--' ,label='analytical solution')


plt.xlabel('function input') # x axis label 
plt.ylabel("f'(x)") #y axis label 
plt.title('forward differencing vs analytical solution') # graph title
plt.legend() # adding a legend with the labels declared in 71-73
plt.show()# showing the graph 

plt.plot(a,centdifs, label='central differencing (.15)') #82-84 plotting the central differencing derivatives and the analytical solution  
plt.plot(a,centdifb, label = "central differencing (.5)") 
plt.plot(a,first, '--', label='analytical solution')


plt.xlabel('function input') # x axis label 
plt.ylabel("f'(x)") #y axis label 
plt.title('central differencing vs analytical solution') # graph title
plt.legend() # adding a legend with the labels declared in 64-68
plt.show()# showing the graph 

#%%
'''
#QUESTION 1
#PART C 
'''


z = np.logspace(-16, 0, 17) # creating a numpy.ndarray of the various values for h 

em = np.finfo(np.float64).eps # computer round off error 

#fde = np.array(abs(ford(1, z)- fd(1)))# forward differencinf error from my calculations
cde = np.array(abs(centd(1, z) - fd(1))) # central differencing error from my calculations
fde = np.array(abs(ford(1, z) - fd(1))) # forward differencing error from my calculations
Efd = np.array(abs((z/2)*(sd(1)) + 2*abs(f(1))*em/z)) #forward differencing error from given error formula 5
Ecd = np.array(abs(((z**2)/24)*(td(1)) + 2*abs(f(1))*em/z)) #forward differencing error from given error formula 13


plt.plot(z, fde, label = 'ε - fd my calculations') # libnes 111-114 : plotting values from lines 105 to 108
plt.plot(z, cde, label = 'ε - cd my calculations')
plt.plot(z, Efd, label ='ε - fd')
plt.plot(z, Ecd, label = 'ε - cd')

plt.xscale('log') # setting the x axis to log scale because h is log values 
plt.yscale('log')  # setting the y axis to log scale becuase the error values are small


plt.xlabel('h value ') # x axis label 
plt.ylabel(" absolute error") #y axis label 
plt.title('Calculated eroor as a function of numerous h values') # graph title
plt.legend() # adding a legend with the labels declared in 111 - 114
plt.show()# showing the graph 


#%%
'''
#QUESTION 1
#PART D
'''

Rfd = np.array(abs(((4*f(1+(z/2)) - f(1+z) - 3*f(1))/z)-fd(1))) # calculating richardson forward differencing error
Rcd = np.array(abs(((8*f(1+(z/4)) + f(1-(z/2)) - f(1+(z/2)) - 8*f(1-(z/4)))/(3*z))-fd(1))) # calculating richardson central differencing error

plt.plot(z, Rfd, label = 'ε - Rich fd') # lines 136 to 139 : plotting values of richardson extrapolation and standard fd, cd formulas  
plt.plot(z, Rcd, label = 'ε - Rich cd') 
plt.plot(z, Efd, label ='ε - fd') 
plt.plot(z, Ecd, label = 'ε - cd')

plt.xscale('log') # setting the x axis to log scale because h is log values 
plt.yscale('log') # setting the y axis to log scale becuase the error values are small


plt.xlabel('h value ') # x axis label 
plt.ylabel(" absolute error") #y axis label 
plt.title('Calculated eroor as a function of numerous h values') # graph title
plt.legend() # adding a legend with the labels declared in lines 136- 139
plt.show()# showing the graph 


#%%
'''
#QUESTION 2
#PART A
'''

# lines 159 to 173 were given 
def legendre(n,x):
    if n==0: # P
        val2 = 1. # P0
        dval2 = 0. # P0'
    elif n==1: # derivatives
        val2 = x # P1'
        dval2 = 1. # P1'
    else:
        val0 = 1.; val1 = x # sep P0 and P2 to start recurrence relation
        for j in range(1,n):
            val2 = ((2*j+1)*x*val1 - j*val0)/(j+1) 
            # P_j+1=(2j+1)xP_j(x)-jP_j-1(x)  / (j+1), starts from P2
            val0, val1 = val1, val2
        dval2 = n*(val0-x*val1)/(1.-x**2) # derivative
    return val2, dval2

def f(n,x) : return (x**2 -1)**n # defining our new function  

h = 0.01 # setting h 

# lines 180 to 194 were only slightly modified from the given code 
def cd_1(f,n,x,h): # defining central differencing for one derivative 
    cd = (f(n,x+h/2) - f(n,x-h/2))/h # calculating the cd
    return cd # definition output 
#cd_2, cd_3 and cd_4 work the same as cd_1 however they reference earlier central differences 
def cd_2 (f ,n ,x , h ):
    cd = ( cd_1 (f ,n , x+h/2 , h ) - cd_1 (f ,n , x-h/2 , h ) )/h
    return cd

def cd_3 (f ,n ,x , h ):
    cd = ( cd_2 (f ,n , x+h/2 , h ) - cd_2 (f ,n , x-h/2 , h ) )/h
    return cd

def cd_4(f, n, x, h):
    cd = ( cd_3(f ,n , x+h/2 , h ) - cd_3(f ,n , x-h/2 , h ) )/h
    return cd


def rodrigues(n , x , h): # putting the central differences from above into a loop for ease of access when plotting 
    if n == 1 :
        return (1/2**n*ma.factorial(n))*(cd_1(f,n,x,h))
    if n == 2 :
        return (1/((2**n)*ma.factorial(n)))*cd_2(f,n,x,h)
    if n == 3 :
        return (1/((2**n)*ma.factorial(n)))*cd_3(f,n,x,h)
    if n ==  4:
      return (1/((2**n)*ma.factorial(n)))*cd_4(f,n,x,h)


# lines 209 to 223 were only slightly modified from the given code
def plotcomp(nsteps):

    xs = [i/nsteps for i in range (-nsteps+1,nsteps)]
    for n in range(1,5): # this will compute P1 to P4
        ys = [legendre(n,x)[0] for x in xs]
        plt.plot(xs, ys, 'k-', label='recurrence n={0}'.format(n), linewidth=3)
        plt.xlabel('$x$', fontsize=20)
        plt.ylabel("$P_n(x)$", fontsize=20)
        ys = [rodrigues(n, x, h) for x in xs] # 217/218 plotting my function 
        plt.plot(xs, ys, 'r--', label= 'rodrigues n = {0}'.format(n), linewidth = 3)
        plt.legend(loc="best")
        plt.show() # create a new graph for each n

nsteps = 200 # number of x points (200 should be enough)
plotcomp(nsteps)


#%%
'''
QUESTION 2
PART B
'''

def new_rodrigues(n,x,h,i): #defining a more efficient rodrigues formula 
    if (i== 1):
        j = (f(n,x+h/2) - f(n,x-h/2))/h
    else : # if i is not 1 then we must continue to take derivatives until i is 1, we do this looping the function through itself and and using output 1 (which is jsut j)
        j = (new_rodrigues(n,x +h/2, h, i-1)[1] - new_rodrigues(n, x-h/2, h, i-1)[1])/h
    return j*(1/((2**n)*ma.factorial(n))), j # output zero is plotted, output 1 is used in the functions loop 

# lines 240 to 254 are only slightly modified from given code
def plotcomp(nsteps):
   xs = [i/nsteps for i in range (-nsteps+1,nsteps)]
   for n in range(1,9): # this will compute P1 to P8
        ys = [legendre(n,x)[0] for x in xs]
        plt.plot(xs, ys, 'k-', label='recurrence n={0}'.format(n), linewidth=3)
        plt.xlabel('$x$', fontsize=20)
        plt.ylabel("$P_n(x)$", fontsize=20)
        ys = [new_rodrigues(n, x, h, n)[0] for x in xs] # 247/248 plotting the more efficient rodrigues formula
        plt.plot(xs, ys, 'r--', label= 'new rodrigues n = {0}'.format(n), linewidth = 3)
        plt.legend(loc="best")
        plt.show() # create a new graph for each n

nsteps = 200 # number of x points (200 should be enough)
plotcomp(nsteps)
