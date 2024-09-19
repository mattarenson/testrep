import numpy as np

def compute_order(x, xstar):
    x = x[:-1]
    diff1 = np.abs(x[1::]-xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()), 1)
    print('the order equation is')
    print('log(|p_{n+1}-p|) = log(lambda) + alpha*log(|p_n-p|) where')
    print('lambda = ' + str(np.exp(fit[1])))
    print('alpha = ' + str(fit[0]))
    return [fit, diff1, diff2]

def driver():
    # Test functions 
    f1 = lambda x: (10/(x + 4))**0.5  # Corrected the power operation

    Nmax = 100
    tol = 1e-10  # Changed tolerance to 10^-10

    # Test f1
    x0 = 1.5
    xstar, ier, x = fixedpt(f1, x0, tol, Nmax)
    print('Fixed Point (f1):', xstar)
    print('f1(xstar):', f1(xstar))
    print('Error message reads:', ier)
    print('Approximations:\n', x) 
    
    order_results = compute_order(x, xstar)
    print('Order of convergence results:', order_results)

# Define routines
def fixedpt(f, x0, tol, Nmax):
    x = []  
    count = 0

    while count < Nmax:
        x1 = f(x0)
        x.append(x1)  #
        
        if abs(x1 - x0) < tol:
            xstar = x1
            ier = 0
            return xstar, ier, np.array(x) 
        
        x0 = x1  
        count += 1 

    xstar = x1
    ier = 1
    return xstar, ier, np.array(x)  # Return as a column vector

driver()


