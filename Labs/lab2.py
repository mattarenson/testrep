import numpy as np

def driver():
    # use routines    
    f = lambda x: np.sin(x)
    a = .5
    b = (3*np.pi)/4
    tol = 1e-5

    [astar, ier] = bisection(f, a, b, tol)
    print('The approximate root is', astar)
    print('The error message reads:', ier)
    print('f(astar) =', f(astar))

# define routines
def bisection(f, a, b, tol):
    
    fa = f(a)
    fb = f(b)
    
    if (fa * fb > 0):
        ier = 1
        astar = a
        return [astar, ier]

    # Verify end points are not a root 
    if (fa == 0):
        astar = a
        ier = 0
        return [astar, ier]

    if (fb == 0):
        astar = b
        ier = 0
        return [astar, ier]

    count = 0
    d = 0.5 * (a + b)
    
    while (abs(d - a) > tol):
        fd = f(d)
        
        if (fd == 0):
            astar = d
            ier = 0
            return [astar, ier]
        
        if (fa * fd < 0):
            b = d
        else: 
            a = d
            fa = fd
        
        d = 0.5 * (a + b)
        count = count + 1
      
    astar = d
    ier = 0
    return [astar, ier]

# Call the driver function
driver()

## Question 1
# a, c return the same ouput, where the root is one, part b doesn't work because the root 1 is outside of the interval. The reason
# 0 is not a root is because the function doesn't cross the x-axis at that root.

## Question 2
# Part a worked as desired with proper accuracy. Part b didnt work due to the function not crossing the x-axis at the root inside the bounds
# For part c the desired root is the on one of the bounds, therefore the accuracy is not as desired as wanted

# import libraries
import numpy as np
    
def driver():

# test functions 
     f1 = lambda x: x*(1+((7-x**5)/(x**2)))
# fixed point is alpha1 = 1.4987....

     f2 = lambda x: x-((x**5-7)/(x**2))
#fixed point is alpha2 = 3.09... 

     Nmax = 100
     tol = 1e-10

# test f1 '''
     x0 = 1
     [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
    
#test f2 '''
     x0 = 1
     [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f2(xstar):',f2(xstar))
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

## Question 3
# Yes 7^(1/5) is a fixed point for all of these functions
# Fixed point doesn't converge for functions a and b and it converges for functions c and d. For 
# a, values far from fixed pt lead to very large numbers, a similar effect is seen in b as the second
# term of the function does the same. However, for c, second term is behaving much less agressively as it
# is in essence 1/x, d's fixed denominator ensures corrections do not become too large