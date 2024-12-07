import numpy as np

## Question 1a
A = np.array([
    [6, 2, 2],
    [2, 2/3, 1/3],
    [1, 2, -1]
])
b = np.array([-2, 1, 0])
c = np.array([2.6, -3.8, -5])

result = A @ c
print(result)

## Question 1b
A = np.array([
    [6, 2, 2, -2],
    [2, 2/3, 1/3, 1],
    [1, 2, -1, 0]
], dtype=float)

def round(x):
    return np.round(x, 4)
n = len(A)
for i in range(n):
    for j in range(i+1, n):
        mult = round(A[j, i] / A[i, i])
        A[j, i:] = round(A[j, i:] - mult * A[i, i:])

solution = np.zeros(n)
for i in range(n-1, -1, -1):
    solution[i] = round((A[i, -1] - np.dot(A[i, i+1:n], solution[i+1:])) / A[i, i])

print(solution)

## Question 1c
A = np.array([
    [6, 2, 2, -2],
    [2, 2/3, 1/3, 1],
    [1, 2, -1, 0]
], dtype=float)

def round(x):
    return np.round(x, 4)

n = len(A)
for i in range(n):
    max= np.argmax(abs(A[i:, i])) + i
    if max != i:
        A[[i, max]] = A[[max, i]]
    for j in range(i+1, n):
        mult = round(A[j, i] / A[i, i])
        A[j, i:] = round(A[j, i:] - mult * A[i, i:])

solution = np.zeros(n)
for i in range(n-1, -1, -1):
    solution[i] = round((A[i, -1] - np.dot(A[i, i+1:n], solution[i+1:])) / A[i, i])

print(solution)

## Question 2

A = np.array([
    [12, 10, 4],
    [10, 8, -5],
    [4, -5, 3]
], dtype=float)

def householder(A):
    x = A[1:, 0]
    norm = np.linalg.norm(x)
    sign = 1 if x[0] >= 0 else -1
    v = x + sign * norm * np.array([1, 0])
    v = v / np.linalg.norm(v)
    H = np.eye(3)
    H[1:, 1:] -= 2 * np.outer(v, v)
    A1 = H @ A @ H.T
    return H, A1

H, A1 = householder(A)

print(A1)

## Question 3a

import time

def power_test():
    # Test for n = 4, 8, 12, 16, 20
    for n in range(4, 21, 4):
        # Generate the Hilbert matrix
        A = np.array([[1 / (i + j - 1) for j in range(1, n+1)] for i in range(1, n+1)])
        
        # Power Method parameters
        Nmax = 100
        start = time.time()
        lam, v, la = powerIt(A, Nmax)
        end = time.time()

        # Print results
        print(f"Eigenvalue: {lam:.6f}")
        print(f"Eigenvector {v[:10]}") 
        print(f"Iterations: {len(la)}")

def powerIt(A, Nmax):
    n = A.shape[0]
    q = np.random.rand(n)
    q = q / np.linalg.norm(q) 
    la = [] 
    for j in range(Nmax):
        z = A @ q
        q = z / np.linalg.norm(z)
        l = q.T @ A @ q  
        la.append(l)
    return l, q, np.array(la)

power_test()

## Question 3b

def power_test():
    n = 16
    A = np.array([[1 / (i + j - 1) for j in range(1, n+1)] for i in range(1, n+1)])
    exact_smallest_eigenvalue = min(np.linalg.eigvals(A))
    A_inv = np.linalg.inv(A)
    Nmax = 100
    start = time.time()
    lam_inv, v, la = powerIt(A_inv, Nmax)
    end = time.time()
    smallest_eigenvalue = 1 / lam_inv
    error = abs(smallest_eigenvalue - exact_smallest_eigenvalue)

    print(f"Smallest Eigenvalue  {smallest_eigenvalue}")
    print(f"Exact Smallest Eigenvalue {exact_smallest_eigenvalue}")

def powerIt(A, Nmax):
    n = A.shape[0]
    q = np.random.rand(n)
    q = q / np.linalg.norm(q)
    la = [] 
    for j in range(Nmax):
        z = A @ q
        q = z / np.linalg.norm(z)  
        l = q.T @ A @ q 
        la.append(l)
    return l, q, np.array(la)

power_test()

## Question 3c

def powerIt(A, Nmax):
    n = A.shape[0]
    q = np.random.rand(n)
    q = q / np.linalg.norm(q)  
    la = []  
    for j in range(Nmax):
        z = A @ q
        q = z / np.linalg.norm(z)  
        l = q.T @ A @ q 
        la.append(l)
    return l, q, np.array(la)

n = 16
A = np.array([[1 / (i + j - 1) for j in range(1, n+1)] for i in range(1, n+1)])

exact_eigenvalues = np.linalg.eigvals(A)
exact_smallest_eigenvalue = min(np.abs(exact_eigenvalues))

A_inv = np.linalg.inv(A)
Nmax = 100
lam_inv, v, _ = powerIt(A_inv, Nmax)
approx_smallest_eigenvalue = 1 / lam_inv

A_comp = np.linalg.inv(A_inv)  
E = A_comp - A 
E_norm = np.linalg.norm(E, 2)

cond_number = np.linalg.norm(A, 2) * np.linalg.norm(A_inv, 2)
error = abs(approx_smallest_eigenvalue - exact_smallest_eigenvalue)
bf_bound = cond_number * E_norm

print(f"Bound: {bf_bound}")


