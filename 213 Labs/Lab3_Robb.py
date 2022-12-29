#%%
'''
Ashley Robb 
Student id: 20203465
ENPH 213 ASSIGNMENT 3
'''
print('\nAshley Robb \n\nStudent id: 20203465 \n\nENPH 213 ASSIGNMENT 3')
#%%
'''
QUESTION 1 
PART A
'''
import numpy as np
#import sympy as sm  
#from sympy import *  
#import cmath as cm
#import numpy as np
#import matplotlib.pyplot as plt
import math as ma

print('\n\n------------- RESULTS FROM QUESTION 1 PART A: -------------\n\n')


def array(n): # defining a function to create an n by n array that starts at 21
    A = np.ones((n,n)) # filling array with ones
    for j in range(0, n): #for loop for column positions 
        h = n*j +21 # element value starting at 21
        for i in range(0,n): #for loop for row position 
            A[j,i] = i + h # summing the variables for element output and position 
    return A # returning the matix

A = array(4) # Declaring A for use in future questions 

print('check the function output:\n', A,'\n') # checking the function output 

def lowtri(mat): # defining a function for the lower triangle method
    matlow = np.zeros((len(mat), len(mat))) # creating an array of zeros the same size as the matrix input 
    for i in range(0, len(mat)): # creating a loop to run through the diaganol function
        matlow[i:,i] = mat [i:,i] # replacing all elements in the diaganol and to the left with the original matrix values 
    return matlow # returning the output 

print('lower triangle of function:\n', lowtri(A), '\n') # printing the new matrix 

def uptri(mat): # defining a function for the upper triangle matrix 
    matup = np.zeros((len(mat), len(mat))) # creating a matrix of zeros 
    for i in range(0, len(mat)): # creating a loop that runs through the diaganol of the function 
        matup[i,i:] = mat [i,i:] # replacing values in the diaganol and to the right with the same elements from the original matrix 
    return matup # returning the new matrix 

print('upper triangle of function:\n', uptri(A), '\n') # printing the new matrix


frob = lambda mat : (ma.sqrt(np.sum(mat**2))) # function for Euclidean norm (or Frobenius norm)

def infin(mat): # defining a function for infinity norm 
    vals = [] # defining an empty list 
    for i in range (0, len(mat)): # creating a loop that repeats for the length of the inoutted matrix 
        vals.append(np.sum(mat[i])) # appending the vals variable to sum the current value of the matrix  
    return max(vals) # returning the maximum element from the vals list

print('\nEuclidean norm (Frobenius norm) from my function: ' , frob(A)) # printing my value from the frob function 

print("\nEuclidean norm (Frobenius norm) from Python's function: " , np.linalg.norm(A)) # printing Python's value from their frob function (same output as mine)

print('\nInfinity norm (maximum row-sum norm) from my function: ', infin(A)) # printing my value from my infin function 

print("\nInfinity norm (maximum row-sum norm) from Python's function:" , np.linalg.norm(A, np.inf)) # printing pythons value form their infinity norm function (same output as mine)

#%%
'''
QUESTION 1
PART B
'''

print('\n\n\n------------- RESULTS FROM QUESTION 1 PART B: -------------\n\n')

# I created two functions for Q1. b) i) one uses lots of numpy documentation as the question asks, however the other is much more logical and more efficient 
def inv_numpy(n): # defining a function to create the given matrix with zeros, ones, and negative ones 
    mat_inv = (-1)*np.ones((n, n)) # creating an array of negative ones the size of the input n 
    mat_ones = np.ones((n,n)) # creating an array of ones the size of the input n 
    mat_zer = np.zeros((n,n)) # creating an array of zeros the size of the input n 
    for i in range(0, n): # cylcing through the array n times 
        mat_inv[i:,i] = mat_zer[i:,i] # changing the array of negative ones to be zeros for the bottom left side 
        mat_inv[i,i] = mat_ones[i,i] # changing the array to be ones for the diagonal 
    return mat_inv

A = inv_numpy(4) # setting matrix A to be the matrix defined above with 4 x 4 size 

print('Matrix A (4 x 4 version) from question 1.b) using defined function with numpy ducumentation: \n', A) # printing the matrix 


def inv(n): # defining a more efficietn function to do the same as the other one expcept it uses less of the numpy documentation 
    mat_inv = np.zeros((n, n))
    for i in range(0, n):
        mat_inv[i,i:] = -1
        mat_inv[i,i] = 1
    return mat_inv

# my second function is very similar to my fisrt, except instead of declaring three array and merging them 
# accordingly, I created one and adjusted it where necesary. This was a shorter function but the question 
# specifically asked to use numpy documentation to create the matrix so I wasn't sure which of my methods 
# was correct for this questions 

print('More efficient/simple function for matrix A (4 x 4 version) from question 1.b) however this one only uses one form of numpy documentation:\n', inv(4), '\n') # printing the second function with size 4 x 4 


def alternating_ones(n): # creating the b array to be alternating ones 
    b = np.ones((n,1)) # setting the array to be of size n but only one column 
    b[1::2, ::2] = -1 # setting every other value to be -1 
    return b

b = alternating_ones(4) # setting b to use the previous matrix with size 4  

print("\nMatrix b of alternating +1's and -1's : \n" , b) # printing b 

print('\nFirst three solutions of matrix x from system of linear equations Ax = b (no perturbation) :\n', np.linalg.solve(A,b)[:3]) # printing the first three elements of x 

def perturb_bl(mat): # definign a function to perturb the bottom left corner of a matrix 
    mat[(len(mat)-1), 0] = mat[(len(mat)-1),0] - 0.001 # subtracting 0.001 from the bottom left corner 
    return mat # returnign the matrix 

A_p = perturb_bl(A) # setting A_p (A perturbed) to be the A matrix with the perturbed corner from the above function 

print('\nMatrix A with perturbed bottom left corner : \n' , A_p) # printing the new matrix A with bottom left erturbation 


print('\nFirst three solutions of matrix x from system of linear equations Ax = b (bottom left perturbation) :\n', np.linalg.solve(A_p,b)[:3]) # printing the first three elements of x 


A = inv_numpy(16) # creating new matrix A of size 16 x 16 

print('Matrix A (16 x 16 version) no puturbation : \n', A) # printign the new matrix 

b = alternating_ones(16) # creating a new array of alternating ones with 16 rows and 1 column 

print("\nMatrix b of alternating +1's and -1's from question 1.b) : \n" , b) # printing the new b matrix 

print('\nFirst three elements of matrix x from system of linear equations Ax = b (n = 16) no purturbation :\n', np.linalg.solve(A,b)[:3]) # printing the first three elements of x 

A_p = perturb_bl(A) # giving the new a value a perturbed corner 

print('\nNew matrix A with perturbed bottom left corner : \n' , A_p) # printing the new A matrix of size 16 x 16 with the perturbed corner 

print('\nFirst three elements of matrix x from system of linear equations Ax = b (with perturbation of -0.001 in the bottom left) :\n', np.linalg.solve(A_p,b)[:3]) # printing the first three elements of x 


#%%
'''
QUESTION 2
PART A
'''

print('\n\n\n------------- RESULTS FROM QUESTION 2 PART A: -------------\n\n')

def backsub1(U, bs): # given backsub function 
    n = bs.size
    xs = np.zeros(n)
    xs[n-1] = bs[n-1]/U[n-1, n-1]
    for i in range(n-2, -1,-1):
        bb = 0
        for j in range (i+1, n):
            bb += U[i, j]*xs[j]
            xs[i] = (bs[i] - bb)/U[i, i]
    return xs

def backsub2(U,b): # my more efficient backsub function 
    n = len(U) # setting n to be the length of the matrix (used for looping purposes)
    x = np.zeros(len(b)) # setting an array x to be the size of inputted array b 
    for i in reversed(range (0, n)): # running through the values of U 
        sums = sum(U[i, (i+1):] * x[(i+1):]) # summing the value of of U times x 
        x[i] = (b[i] - sums) / U[i,i] # updating x 
    return x # returning an joutput 

U = uptri(array(5000)) # setting U to be an upper triangle array that is 5000 x 5000 staring at a value of 21 
b = U[0] # setting b to be the first row of U 

from timeit import default_timer # timer function given 
def time(function, a , b):
    timer_start = default_timer ()
    function(a,b)
    timer_end = default_timer ()
    time1 = timer_end - timer_start
    return time1
    
print('First three elements of given/inneficient back sub : \n' , backsub1(U, b)[:3] , '\nand its run time : ', time(backsub1, U,b))
# printing the first three values and time of the given backsub function 

print('\n\nFirst three elements of my more efficient back sub :\n' , backsub2(U,b)[:3] , '\nand its run time :' , time(backsub2, U,b))
#printing the first three values and the time of my backsub function 

#%%
'''
QUESTION 2
PART B
'''

print('\n\n\n------------- RESULTS FROM QUESTION 2 PART B: -------------\n\n')

def gaus(mat , b): # writing a function for Gaussian Elmination 
    matupdiag = np.copy(mat) # copying the inputted matrix  
    bgaus = np.copy(b) # copying the inputted b mastrix 
    for i in range(0, len(mat)): # for loop to run through the length of the inputed matrix 
        for j in range(0,i): # nested for loop to run through the length of the inputted matrix again 
            coef = matupdiag[i, j] / matupdiag[j, j] # defining the coefficient for the guassian process 
            matupdiag[i] = matupdiag[i] - coef*(matupdiag[j]) #  updating the matrix 
            bgaus[i] = bgaus[i] - coef*bgaus[j] # updating the b value 
    Gaussian = backsub2(matupdiag , bgaus) # implementing my backsub function 
    return Gaussian # outputting my result 

A = np.ones((3,3)) # creating a three by three array of ones
A[0,0] = 2; A[1,2] = -2; A[2,1] = 2 # changing values in the array to match the given array 
b =np.array([8,-2,2]) # creating the b array that is given 

print('The Gaussian Elimination solution x to the given A and b is: ' , gaus(A,b)) # printing solution of the function 

#%%
'''
QUESTION 2
PART C
'''

print('\n\n\n------------- RESULTS FROM QUESTION 2 PART C: -------------\n\n')

A = np.ones((3,3)) # creating a 3 by 3 array of ones 
A[0,0] = A[1,0] = A[2,1] = 2 # changing elements to match the given array 
A[1,2] = -4

B =np.array([8,-2,2]) # creating the given b array 

def gaus_partial_pivot(mat , b): # defining a funciton to use gaussian with partial pivotitng 
    matpp = np.copy(mat) # copying the inputted matrix to a new matrix partial pivoting matrix 
    bpp = np.copy(b) # copying b to a new partial pivoting array 
    for i in range(0, len(mat)): # running a loop to run as many times as the length of the inputted matrix 
        n = len(matpp) - 1 # setting n to be one less than the length of the partial pivot matrix 
        k = np.argmax(np.abs(matpp[:,n - i])) # setting k to be the max value of the column 
        matpp[[0,k]] = matpp[[k,0]] # switching the rows based on k 
        bpp[[0,k]] = bpp[[k , 0]] # switching the elements in b to match the switch in the matrix 
    return(gaus(matpp , bpp))   # calling the gaus function from above and returning it with the new values

print('Results after modifying function to partially pivot :\n' ,gaus_partial_pivot(A , B)) # printing the gausian results of the partially pivotted matrix 


#%%
'''
QUESTION 2
PART D
'''

print('\n\n\n------------- RESULTS FROM QUESTION 2 PART D: -------------\n\n')

A = np.zeros((3,3)) # setting an array of sizes 3 by 3 to have all elements of zero 
A[0 , 0] = 1 ; A[0 , 1] = 2 ; A[0 , 2] = 3 # adjusting the array to match the given array 
A[1 , 1] = 1 ; A[1 , 2] = 4
A[2 , 0] = 5 ; A[2 , 1] = 6 

b = np.zeros((3,3)) # setting an arrat of 3 by 3 to be all zeros 
b[0,0] = b[1,1] = b[2,2] = 1 # adjusting the array to be the identity matrix 


x = gaus_partial_pivot(A , b[:,0]) # setting x to be the output of the gaussian partial pivot with the matrix A and the first column of b 

y = gaus_partial_pivot(A , b[:,1]) # setting y to be the output of the gaussian partial pivot with the matrix A and the second column of b 

z = gaus_partial_pivot(A , b[:,2]) # setting z to be the output of the gaussian partial pivot with the matrix A and the third column of b 

A_inverse = np.zeros((3,3)) # setting A inverse to be an array of 3 by 3 zeros 

for i in range(0, 3): # adjusting the values fo A inverse with a for loop 
    A_inverse[i , 0] = x[i] # setting the first column of A inverse to be x 
    A_inverse[i , 1] = y[i] # setting the second column of A inverse to be y 
    A_inverse[i , 2] = z[i] # setting the third column of A inverse to be z

print('The inverse of matrix A is :\n\n' , A_inverse) # printing the inverse matrix 

eq_8 = np.matmul( A , A_inverse) # using equation 8 to check my solutions 

print('\n\nChecking solution with equation 8 : \n\n' , eq_8) # printing my check from equation 8

print('\n\n***please note: my check is slightly off, I believe that this is due to machine error, if you look closly each value is only of by e-15 or e-16') # explaining my findings 
    




































  