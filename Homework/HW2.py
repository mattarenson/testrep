## Q3
import numpy as np
import math as math
f = lambda x: (math.exp(x))-1

result = f(9.999999995000000*10**(-10))
print(result)

## part d

eps = 1*10**-16
x = 9.999999995*10**-10

n = 0
term = x
while abs(term) >= eps:
    n += 1
    term = ((x)**n)/math.factorial(n)
    
print(n)

## Question 4 part a

increment = np.pi / 30
t = np.arange(0, np.pi + increment, increment)

y = np.cos(t)

def sumof(k, N):
    result = np.sum(t[k-1:N] * y[k-1:N])
    print(f"The sum is: {result}")

k = 1
N = len(t)

sumof(k, N)

# Q4 Part b

import matplotlib.pyplot as plt
import random as random

R = 1.2
delta_r = .1
f = 15
p = 0

theta = np.linspace(0, 2*np.pi, 100)

x_theta = R*(1+delta_r*np.sin(f*theta+p))*np.cos(theta)
y_theta = R*(1+delta_r*np.sin(f*theta+p))*np.sin(theta)

plt.plot(x_theta, y_theta)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()


for i in range(1, 11):
    R = i
    delta_r = 0.05
    f = 2+i
    p = random.uniform(0, 2)
    x_theta = R*(1+delta_r*np.sin(f*theta+p))*np.cos(theta)
    y_theta = R*(1+delta_r*np.sin(f*theta+p))*np.sin(theta)
    plt.plot(x_theta, y_theta)
    plt.gca().set_aspect('equal', adjustable='box')
plt.show()

