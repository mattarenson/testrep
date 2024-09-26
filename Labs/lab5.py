import numpy as np

def driver():
    f = lambda x: x**3 + x - 4
    fp = lambda x: 3*x**2 + 1
    ga = lambda x: (f(x)) * (6*x) / (fp(x))  
    a = 1
    b = 4
    tol = 1e-7
    Nmax = 100  

    [astar, ier] = bisection(f, ga, a, b)
    print('The approximate root is', astar)
    print('The error message reads:', ier)
    print('Bounds =', [astar])

    
    if ier == 0:
        p0 = astar
        [p, pstar, info, it] = newton(f, fp, p0, tol, Nmax)
        print('Newton Method: The approximate root is', pstar)
        print('Newton Method: Iterations:', it)
        print('Newton Method: Error message:', info)

def bisection(f, ga, a, b):
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        ier = 1
        astar = a
        return [astar, ier]

    if fa == 0:
        return [a, 0]
    if fb == 0:
        return [b, 0]

    count = 0
    d = 0.5 * (a + b)
    while abs(ga(d)) > 1:
        fd = f(d)
        if fd == 0:
            return [d, 0]
        if fa * fd < 0:
            b = d
        else:
            a = d
            fa = fd
        d = 0.5 * (a + b)
        count += 1

    astar = d
    ier = 0
    return [astar, ier]

def newton(f, fp, p0, tol, Nmax):
    p = np.zeros(Nmax + 1)
    p[0] = p0
    for it in range(Nmax):
        p1 = p0 - f(p0) / fp(p0)
        p[it + 1] = p1
        if abs(p1 - p0) < tol:
            return [p, p1, 0, it]
        p0 = p1
    return [p, p1, 1, it]

# Run the driver function
driver()

