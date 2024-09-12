import numpy as np
import numpy.linalg as la
import math
def driver():
    n = 2
    y = [1, 0]
    w = [0, 1]

    dp = dotProduct(y,w,n)

    print('the dot product is : ', dp)
    return
def dotProduct(x,y,n):
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp
driver()

def matrixMultiplication(mat, x):
    new = []
    i = len(mat)
    j = len(mat[0])
    for i in range(i):
        sum = 0
        for j in range(j):
            sum += mat[i][j] * x[j]
        new.append(sum)
    return new

    
