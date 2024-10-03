import numpy as np

f = np.cos
s = np.pi / 2

h_values = [0.01 * (2 ** -i) for i in range(10)]

h_values = np.array(h_values)

forward_diff = (f(s + h_values) - f(s)) / h_values
centered_diff = (f(s + h_values) - f(s - h_values)) / (2 * h_values)

print("Forward difference approximations:", forward_diff)
print("Centered difference approximations:", centered_diff)

true_derivative = -np.sin(s)

forward_errors = np.abs(forward_diff - true_derivative)
centered_errors = np.abs(centered_diff - true_derivative)

forward_orders = np.log(forward_errors[1:] / forward_errors[:-1]) / np.log(h_values[1:] / h_values[:-1])
centered_orders = np.log(centered_errors[1:] / centered_errors[:-1]) / np.log(h_values[1:] / h_values[:-1])

print(forward_orders)
print(centered_orders)

## order of forward is 1 and order of centered is 2

import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():

    x0 = np.array([1, 0])
    
    Nmax = 100
    tol = 1e-10

    t = time.time()
    for j in range(20):
      [xstar,ier,its] =  SlackerNewton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Lazy Newton: the error message reads:',ier)
    print('Lazy Newton: took this many seconds:',elapsed/20)
    print('Lazy Newton: number of iterations is:',its)
     
     
def evalF(x): 

    F = np.zeros(2)
    
    F[0] = 4*x[0]+x[1]**2-4
    F[1] = x[0]+x[1]+np.sin(x[0]-x[1])
    return F
    
def evalJ(x): 

    
    J = np.array([[4, 2*x[1]], 
        [1-np.cos(x[0]-x[1]), 1-np.cos(x[0]-x[1])]])
    return J

           
def SlackerNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ(x0)
    Jinv = inv(J)
    xstart = x0
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its]
           
       x0 = x1
       if (np.linalg.norm(x1-xstart)>1):
           J = evalJ(x1)
           Jinv = inv(J)
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]   
        
if __name__ == '__main__':
    # run the drivers only if this is called from the command line
    driver()       

## condition I decided to put is if the difference between x1 and x0 is greater than 1
## Partner got 4 iterations because used .001 as difference I got 19 because I used 1


