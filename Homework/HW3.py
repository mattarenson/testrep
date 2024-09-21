# import libraries
import numpy as np
import matplotlib.pyplot as plt

def driver():
    # Use the expanded version of the function
    f = lambda x: x**3+x-4
    a = 1
    b = 4
    tol = 1e-3

    [astar, ier], count = bisection(f, a, b, tol)
    print('The approximate root is', astar)
    print('The error message reads:', ier)
    print('f(astar) =', f(astar))
    print(count)

def bisection(f, a, b, tol):
    fa = f(a)
    fb = f(b)
    if (fa * fb > 0):
        ier = 1
        astar = a
        return [astar, ier]

    if (fa == 0):
        astar = a
        ier = 0
        return [astar, ier]

    if (fb == 0):
        astar = b
        ier = 0
        return [astar, ier]
    count = 0
    while (b - a) > tol:
        d = 0.5 * (a + b)
        fd = f(d)
        count +=1
        if (fd == 0):
            astar = d
            ier = 0
            return [astar, ier]

        if (fa * fd < 0):
            b = d
        else: 
            a = d
            fa = fd

    astar = 0.5 * (a + b)
    ier = 0
    return [astar, ier], count

driver()


x = np.linspace(-2, 6)
def func1(x):
    return x-4*np.sin(2*x)-3

plt.plot(x, func1(x))
plt.axhline(0, color = 'red', linestyle = '--')
plt.show()

def driver():

# test functions 
     f1 = lambda x: -np.sin(2 * x) + (5 * x / 4) - (3 / 4)



     Nmax = 100
     tol = 1e-10

     x0 = 5
     [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
    



# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    

driver()
