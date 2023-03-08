import numpy as np


def lu_pivot(A):
    n = len(A)
    L = [[0 for i in range(n)] for j in range(n)]
    U = [[0 for i in range(n)] for j in range(n)]
    P = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    for i in range(n):
        pivot_row = max(range(i,n), key=lambda k: abs(A[k][i]))
        if pivot_row != i:
            A[i], A[pivot_row] = A[pivot_row], A[i]
            P[i], P[pivot_row] = P[pivot_row], P[i]
        U[i][i] = A[i][i] - sum(L[i][k] * U[k][i] for k in range(i))
        for j in range(i+1,n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]
        for j in range(i+1,n):
            U[i][j] = (A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))) / L[i][i]
    return [L, U, P]



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
A = np.array([[1, 0, 0, 1, 0],
              [0, 1, 1, 0, 0],
              [0, 0, 1, -1, -1],
              [2, 0, 0, -8, 10],
              [0, 4, -6, 0, -10]])


# Initialize vector b
b = np.array([1, 1, 0, 0, 0])

#output part
print('\nThe solution for different I(s) are given below:\n')
for i in range (len(b)):
    print('I{} = {}'.format(i+1,LU_solver(A,b)[i]))


print('\n\nThe decomposited L and U are:\n')
print('L = {}\n\nU = {}'.format(LU_decomp(A)[0],LU_decomp(A)[1]))
          
