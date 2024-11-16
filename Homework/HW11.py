import numpy as np
from scipy.integrate import quad

def trapezoidal(f, a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return h * (0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1])

def simpsons(f, a, b, n):
    if n % 2 != 0:
        return
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)

    return (h/3) * (y[0]+4 * np.sum(y[1:-1:2])+ 2 * np.sum(y[2:-2:2])+ y[-1])

def main():
    a = -5
    b = 5
    ntrap = 42
    nsimp = 12

    f = lambda s: 1 / (1 + s**2)

    trap_estimate = trapezoidal(f, a, b, ntrap)
    simp_estimate = simpsons(f, a, b, nsimp)
    quad_default, _ = quad(f, a, b, epsabs=1e-6, full_output=False)
    quad_tolerance, _ = quad(f, a, b, epsabs=1e-4, full_output=False)
    print(f"Quad Tolerance 1e-6: {quad_default:.6f}")
    print(f"Quad Tolerance 1e-4: {quad_tolerance:.6f}")
    print(f"Trapezoidal: {trap_estimate}")
    print(f"Simpsons: {simp_estimate}")

if __name__ == "__main__":
    main()

## Question 2

import numpy as np
def simpsons(f, a, b, n):
    if n % 2 != 0:
        return
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)

    return (h/3) * (y[0]+4 * np.sum(y[1:-1:2])+ 2 * np.sum(y[2:-2:2])+ y[-1])

def main():
    a = 1e-5
    b = 1
    nsimp = 4

    f = lambda s: np.cos(1/s)/s

    simp_estimate = simpsons(f, a, b, nsimp)
    print(f"Simpsons: {simp_estimate}")

if __name__ == "__main__":
    main()