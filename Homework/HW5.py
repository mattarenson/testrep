import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

# Question 1a

f = lambda x, y: 3*x**2 - y**2
g = lambda x, y: 3*x*y**2 - x**3 - 1

x = 1.0
y = 1.0

tol = 1e-10

A = np.array([[1/6, 1/18],
                  [0, 1/6]])

maxN = 100
for n in range(maxN):
    F = np.array([f(x, y), g(x, y)])

    G = A @ F
    x1 = x - G[0]
    y1 = y - G[1]
    
    if np.linalg.norm([x1 - x, y1- y]) < tol:
        break

    ier = 0
    x, y = x1, y1
else:
    ier = 1

print(f"Error = {ier}")
print(f"x = {x}, y = {y}")
print(f"f(x, y) = {f(x, y)}")
print(f"g(x, y) = {g(x, y)}")
print(f"Number of iteration= {n}")

# Question 1c

def driver():

    x0 = np.array([1, 1])
    
    Nmax = 100
    tol = 1e-10
    
    t = time.time()
    for j in range(50):
      [xstar,ier,its] =  Newton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Newton: the error message reads:',ier) 
    print('Newton: took this many seconds:',elapsed/50)
    print('Netwon: number of iterations is:',its)
     
def evalF(x): 

    F = np.zeros(2)
    
    F[0] = 3*x[0]**2-x[1]**2
    F[1] = 3*x[0]*x[1]**2-x[0]**3-1
    return F
    
def evalJ(x): 
    J = np.array([[6*x[0], -2*x[1]], 
        [3*x[1]**2-3*x[0]**2, 6*x[0]*x[1]]])
    return J


def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

driver()

# Question 3b

f = lambda x, y, z: x**2 + 4*y**2 + 4*z**2 - 16
fx = lambda x: 2*x
fy = lambda y: 8*y
fz = lambda z: 8*z

x = 1.0
y = 1.0
z = 1.0

tol = 1e-10
ier = 0

max_iterations = 100

for n in range(max_iterations):
    F = np.array([f(x, y, z)])

    fx1 = fx(x)
    fy1 = fy(y)
    fz1 = fz(z)
 
    gradient = fx1**2 + fy1**2 + fz1**2

    x1 = x - (F[0] * fx1) / gradient
    y1 = y - (F[0] * fy1) / gradient
    z1 = z - (F[0] * fz1) / gradient
  
    if np.linalg.norm([x1 - x, y1 - y, z1 - z]) < tol:
        break
    
    ier = 0
    x, y, z = x1, y1, z1
else:
    ier = 1

print(f"Error = {ier}")
print(f"x = {x}, y = {y}, z = {z}")
print(f"Number of iteration = {n}")


