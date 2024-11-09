import numpy as np
import matplotlib.pyplot as plt

def pade_3_3_parta(x):
    numerator = x - (7/60) * x**3
    denominator = 1 + (1/20) * x**2
    return numerator / denominator

def maclaurin_polynomial(x):
    return x - x**3 / 6 + x**5 / 120

def sin(x):
    return np.sin(x)

x_values = np.linspace(0, 5, 500)

sin_values = sin(x_values)
maclaurin_values = maclaurin_polynomial(x_values)
pade_3_3_values_case1 = pade_3_3_parta(x_values)

error_maclaurin = np.abs(sin_values - maclaurin_values)
error_pade_3_3_case1 = np.abs(sin_values - pade_3_3_values_case1)

# Plotting errors
plt.figure(figsize=(10, 6))
plt.plot(x_values, error_maclaurin, label="Maclaurin (6th)", linestyle="--")
plt.plot(x_values, error_pade_3_3_case1, label="Padé(3,3) (part a, c3 = 0)")
plt.xlabel("x")
plt.ylabel("Error")
plt.title("Padé(3,3) Numerator and Denominator Cubic and Maclaurin Polynomial for sin(x)")
plt.legend()

from sympy import symbols, series, Rational, simplify, collect, Eq, solve

x = symbols('x')
a0, a1, a2 = symbols('a0 a1 a2')
b1, b2, b3, b4 = symbols('b1 b2 b3 b4')

taylor_series = x - Rational(1, 6)*x**3 + Rational(1, 120)*x**5
numerator_2_4 = a0 + a1 * x + a2 * x**2
denominator_2_4 = 1 + b1 * x + b2 * x**2 + b3 * x**3 + b4 * x**4
pade_2_4_approximation = numerator_2_4 / denominator_2_4

pade_2_4_series = series(pade_2_4_approximation, x, 0, 6).removeO()

equations_2_4 = [collect(simplify(taylor_series - pade_2_4_series), x).coeff(x, i).simplify() for i in range(6)]

solution_2_4 = solve(equations_2_4, (a0, a1, a2, b1, b2, b3, b4))
solution_2_4

x_values = np.linspace(0, 5, 500)

def sin_function(x):
    return np.sin(x)

def maclaurin_polynomial(x):
    return x - x**3 / 6 + x**5 / 120


def pade_2_4(x):
    numerator = x  
    denominator = 1 + (1/6) * x**2 + (7/360) * x**4
    return numerator / denominator

sin_values = sin_function(x_values)
maclaurin_values = maclaurin_polynomial(x_values)
pade_2_4_values_case_b3_zero = pade_2_4(x_values)

error_maclaurin = np.abs(sin_values - maclaurin_values)
error_pade_2_4_case_b3_zero = np.abs(sin_values - pade_2_4_values_case_b3_zero)


plt.plot(x_values, error_maclaurin, label="Maclaurin (6th)", linestyle="--")
plt.plot(x_values, error_pade_2_4_case_b3_zero, label="Padé(2,4) part b with b3 = 0")
plt.xlabel("x")
plt.ylabel("Error")
plt.title("Padé(2,4) and Maclaurin Polynomial for sin(x)")
plt.legend()
plt.show()

def pade_4_2_case_c1_zero(x):
    numerator = x - (7/60) * x**3  
    denominator = 1 + (1/20) * x**2  
    return numerator / denominator

sin_values = sin_function(x_values)
maclaurin_values = maclaurin_polynomial(x_values)
pade_4_2_values_case_c1_zero = pade_4_2_case_c1_zero(x_values)

error_maclaurin = np.abs(sin_values - maclaurin_values)
error_pade_4_2_case_c1_zero = np.abs(sin_values - pade_4_2_values_case_c1_zero)

plt.plot(x_values, error_maclaurin, label="Maclaurin (6th)", linestyle="--")
plt.plot(x_values, error_pade_4_2_case_c1_zero, label="Padé(4,2) part c with c1 = 0")
plt.xlabel("x")
plt.ylabel("Error")
plt.title("Padé(4,2) and Maclaurin Polynomial for sin(x)")
plt.legend()
plt.grid(True)
plt.show()

## Question 2
from sympy import symbols, Eq, solve

x0, x1, c1 = symbols('x0 x1 c1')
eq1 = Eq(1, 1/2 + c1)
eq2 = Eq(1/2, 1/2 * x0 + c1 * x1)
eq3 = Eq(1/3, 1/2 * x0**2 + c1 * x1**2)
eq4 = Eq(1/4, 1/2 * x0**3 + c1 * x1**3)

solution = solve((eq1, eq2, eq3, eq4), (x0, x1, c1))

print(solution)

