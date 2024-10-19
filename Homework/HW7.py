import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

def f(x):
    return 1/(1+(10*x)**2)

def interpolation(N):
    h = 2/(N-1)
    xpts = np.array([-1+(i+1)*h for i in range(1, N+1)])
    ypts = f(xpts)

    V = np.vander(xpts, increasing=True)

    c = np.linalg.solve(V, ypts)

    x2 = np.linspace(-1, 1, 1001)
    f2 = f(x2)
    p2 = np.polyval(c[::-1], x2)

    plt.plot(xpts, ypts, 'o')
    plt.plot(x2, f2, label='f(x)')
    plt.plot(x2, p2, label='Interpolation')
    plt.title(f'Interpolation with N = {N}')
    plt.xlim(-1, 1)
    plt.ylim(0, 1.5)
    plt.legend()
    plt.grid()
    plt.show()

for N in range(2, 21):
    interpolation(N)

## Question 2

def f(x):
    return 1 / (1 + (10 * x) ** 2)

def weights(xpts):
    n = len(xpts)
    w = np.ones(n) 
    for j in range(n):
        w[j] =1 
        for i in range(n):
            if i != j:
                w[j] *= (xpts[j] - xpts[i])
    return 1 / w  

def barycentric(x, xpts, ypts, w):
    n = len(xpts)
    num = np.zeros(len(x))
    denom = np.zeros(len(x))

    for j in range(n):
        cond = (x != xpts[j]) 
        arc = w[j] / (x[cond] - xpts[j]) 
        num[cond] += arc * ypts[j]  
        denom[cond] += arc  

        with np.errstate(divide='ignore', invalid='ignore'):
            val = np.where(denom != 0, num / denom, 0) 

    return val

def interpolation(N):
    h = 2 / (N - 1)
    xpts = np.array([-1 + i * h for i in range(N)]) 
    ypts = f(xpts)

    w = weights(xpts)

    x2 = np.linspace(-1, 1, 1001) 
    f2 = f(x2)
    p2 = barycentric(x2, xpts, ypts, w)

    plt.figure(figsize=(10, 6))
    plt.plot(xpts, ypts, 'o', label='Data Points')
    plt.plot(x2, f2, label='f(x)')
    plt.plot(x2, p2, label='Barycentric Interpolation')
    plt.title(f'Barycentric Interpolation with N = {N}')
    plt.xlim(-1, 1)
    plt.ylim(0, 1.5)
    plt.legend()
    plt.grid()
    plt.show()

for N in range(2, 21):
    interpolation(N)
  
## Question 3 (a)

def f(x):
    return 1/(1+(10*x)**2)

def interpolation(N):
    xpts = np.cos((2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N))
    ypts = f(xpts)

    V = np.vander(xpts, increasing=True)

    c = np.linalg.solve(V, ypts)

    x2 = np.linspace(-1, 1, 1001)
    f2 = f(x2)
    p2 = np.polyval(c[::-1], x2)

    plt.plot(xpts, ypts, 'o')
    plt.plot(x2, f2, label='f(x)')
    plt.plot(x2, p2, label='Interpolation')
    plt.title(f'Interpolation with N = {N}')
    plt.xlim(-1, 1)
    plt.ylim(0, 1.5)
    plt.legend()
    plt.grid()
    plt.show()

for N in range(2, 21):
    interpolation(N)

## Question 3 (b)

def f(x):
    return 1 / (1 + (10 * x) ** 2)

def weights(xpts):
    n = len(xpts)
    w = np.ones(n) 
    for j in range(n):
        w[j] =1 
        for i in range(n):
            if i != j:
                w[j] *= (xpts[j] - xpts[i])
    return 1 / w  

def barycentric(x, xpts, ypts, w):
    n = len(xpts)
    num = np.zeros(len(x))
    denom = np.zeros(len(x))

    for j in range(n):
        cond = (x != xpts[j]) 
        arc = w[j] / (x[cond] - xpts[j]) 
        num[cond] += arc * ypts[j]  
        denom[cond] += arc  

        with np.errstate(divide='ignore', invalid='ignore'):
            val = np.where(denom != 0, num / denom, 0) 

    return val

def interpolation(N):
    h = 2 / (N - 1)
    xpts = np.cos((2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N))
    ypts = f(xpts)

    w = weights(xpts)

    x2 = np.linspace(-1, 1, 1001) 
    f2 = f(x2)
    p2 = barycentric(x2, xpts, ypts, w)

    plt.figure(figsize=(10, 6))
    plt.plot(xpts, ypts, 'o', label='Data Points')
    plt.plot(x2, f2, label='f(x)')
    plt.plot(x2, p2, label='Barycentric Interpolation')
    plt.title(f'Barycentric Interpolation with N = {N}')
    plt.xlim(-1, 1)
    plt.ylim(0, 1.5)
    plt.legend()
    plt.grid()
    plt.show()

for N in range(2, 21):
    interpolation(N)