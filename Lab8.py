import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.linalg import inv 

def driver():
    f = lambda x: 1 / (1 + (10 * x) ** 2)
    a = -1
    b = 1
    
    # Create points you want to evaluate at
    Neval = 100
    xeval = np.linspace(a, b, Neval)
    
    # Number of intervals
    Nint = 10
    
    # Evaluate the linear spline
    yeval = eval_lin_spline(xeval, a, b, f, Nint)
    
    # Evaluate f at the evaluation points
    fex = np.zeros(Neval)
    for j in range(Neval):
        fex[j] = f(xeval[j])  # Use square brackets for indexing
    
    # Plotting
    plt.figure()
    plt.plot(xeval, fex, label='f(x)', color='blue')
    plt.plot(xeval, yeval, label='Linear Spline', color='orange')
    plt.legend()
    plt.show()
    
    # Calculate error
    err = abs(yeval - fex)
    plt.figure()
    plt.plot(xeval, err, label='Error', color='red')
    plt.legend()
    plt.show()

def linear(x0, y0, x1, y1, alpha):
    m = (y1 - y0) / (x1 - x0)
    b = y0 - m * x0
    y = m * alpha + b
    return y

def eval_lin_spline(xeval, a, b, f, Nint):
    # Create the intervals for piecewise approximations
    xint = np.linspace(a, b, Nint + 1)
    
    # Create vector to store the evaluation of the linear splines
    yeval = np.zeros(len(xeval))
    
    for jint in range(Nint):
        # Find indices of xeval in interval (xint[jint], xint[jint + 1])
        ind = np.where((xeval >= xint[jint]) & (xeval <= xint[jint + 1]))[0]
        
        # Temporarily store your info for creating a line in the interval of interest
        a1 = xint[jint]
        fa1 = f(a1)
        b1 = xint[jint + 1]
        fb1 = f(b1)
        
        # Evaluate the linear spline for the found indices
        for kk in ind:
            yeval[kk] = linear(a1, fa1, b1, fb1, xeval[kk])
    
    return yeval

if __name__ == '__main__':
    # Run the driver only if this is called from the command line
    driver()

def driver():
    
    f = lambda x: np.exp(x)
    a = 0
    b = 1
    
    ''' number of intervals'''
    Nint = 3
    xint = np.linspace(a,b,Nint+1)
    yint = f(xint)

    ''' create points you want to evaluate at'''
    Neval = 100
    xeval =  np.linspace(xint[0],xint[Nint],Neval+1)

#   Create the coefficients for the natural spline    
    (M,C,D) = create_natural_spline(yint,xint,Nint)

#  evaluate the cubic spline     
    yeval = eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D)
    
    
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
        
    nerr = norm(fex-yeval)
    print('nerr = ', nerr)
    
    plt.figure()    
    plt.plot(xeval,fex,'ro-',label='exact function')
    plt.plot(xeval,yeval,'bs--',label='natural spline') 
    plt.legend
    plt.show()
     
    err = abs(yeval-fex)
    plt.figure() 
    plt.semilogy(xeval,err,'ro--',label='absolute error')
    plt.legend()
    plt.show()
    
def create_natural_spline(yint,xint,N):

#    create the right  hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)
    h[0] = xint[i]-xint[i-1]  
    for i in range(1,N):
       h[i] = xint[i+1] - xint[i]
       b[i] = (yint[i+1]-yint[i])/h[i] - (yint[i]-yint[i-1])/h[i-1]

#  create the matrix A so you can solve for the M values
    A = np.zeros((N+1,N+1))
    for i in range(N+1):
        for j in range(N+1):
            if i == j:
                A[i][j] = 4
            if i == j+1 or i + 1 ==j:
                A[i][j] = 1
    
#  Invert A    
    Ainv = inv.A

# solver for M    
    M  = np.linalg.solve
    
#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = # find the C coefficients
       D[j] = # find the D coefficients
    return(M,C,D)
       
def eval_local_spline(xeval,xi,xip,Mi,Mip,C,D):
# Evaluates the local spline as defined in class
# xip = x_{i+1}; xi = x_i
# Mip = M_{i+1}; Mi = M_i

    hi = xip-xi
   
    yeval = 
    return yeval 
    
    
def  eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    
    yeval = np.zeros(Neval+1)
    
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        
#   find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]

# evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
#   copy into yeval
        yeval[ind] = yloc

    return(yeval)
           
driver()             
