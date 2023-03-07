import numpy as np

def LU_decomp(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    # Decomposition
    for i in range(n):
        for j in range(i, n):
            sum = 0
            for k in range(i):
                sum += L[i][k] * U[k][j]
            U[i][j] = A[i][j] - sum

        for j in range(i, n):
            if i == j:
                L[i][i] = 1
            else:
                sum = 0
                for k in range(i):
                    sum += L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - sum) / U[i][i]

    return [L, U]

def LU_solver(A, b):
    L = LU_decomp(A)[0]
    U = LU_decomp(A)[1]
    n = len(A)
    y = np.zeros(n)
    x = np.zeros(n)

    # Solving Ly=b
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += L[i][j] * y[j]
        y[i] = b[i] - sum

    # Solving Ux=y
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum += U[i][j] * x[j]
        x[i] = (y[i] - sum) / U[i][i]

    return x


#Real problem
A = np.array([[1, 0, 0, 1, 0, 0,0],
              [0, 1, 1, 0, 0, 0,0],
              [0, 0, 1, -1, -1, 0,0],
              [1, -1, 0, 0, -1, 0,0],
              [2, 0, 0, -8, 10, 0,0],
              [0, 4, -6, 0, -10, 0,0],
              [2, 4, -6, -8, 0,0,0]])

# Initialize vector b
b = np.array([1, 1, 0, 0, 0, 0, 0])

print(LU_solver(A,b))