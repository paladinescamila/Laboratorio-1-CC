# UNIDAD 2: SISTEMAS DE ECUACIONES LINEALES

import numpy as np
import time

# Método de Sustitución Sucesiva hacia atrás
def sucesiva_hacia_atras(A, b):
    """
    Entrada: una matriz triangular superior A y un vector b.
    Salida: un vector x tal que Ax = b.
    """
    n = len(b) - 1
    x = [None for _ in range(n)] + [b[n] / A[n][n]]
    for i in range(n, -1, -1):
        sumatoria = 0
        for j in range(i+1, n+1):
            sumatoria += A[i][j] * x[j]
        x[i] = (b[i] - sumatoria) / A[i][i]

    return x


# Método de Sustitución Sucesiva hacia adelante
def sucesiva_hacia_adelante(A, b):
    """
    Entrada: una matriz triangular inferior A y un vector b.
    Salida: un vector x tal que Ax = b.
    """
    n = len(b) - 1
    x = [b[0] / A[0][0]] + [None for _ in range(n)]
    for i in range(1, n+1):
        sumatoria = 0
        for j in range(i):
            sumatoria += A[i][j] * x[j]
        x[i] = (b[i] - sumatoria) / A[i][i]

    return x


# Construye la matriz de eliminación para una columna
def matriz_de_eliminacion(A, k, g):
    """
    Entrada: una matriz cuadrada A, un entero k y un booleano g.
    Salida: matriz de eliminación de Gauss (si g es verdadero) o matriz de 
            eliminación de Gauss-Jordan (si g es falso) para la columna Ak.
    """
    n = len(A)
    M = np.identity(n)
    for i in range(k+1, n):
        M[i][k] = (-1) * A[i][k] / A[k][k]
    if (not g):
        for i in range(k):
            M[i][k] = (-1) * A[i][k] / A[k][k]
    
    return M


def row_swap_mat(i, j, n):
    P = np.eye(n)
    P[i] = 0
    P[j] = 0
    P[i][j] = 1
    P[j][i] = 1
    return P

def permutar(A, b, k):
    n = len(A)
    if (k != n-1):
        i = k+1
        while (i != n-1 and A[k][k] == 0):
            A = row_swap_mat(k,i,n).dot(A)
            b = row_swap_mat(k,i,n).dot(b)
            i += 1
            k += 1
    return A, b


# Método de Eliminación de Gauss
def gauss(A, b):
    """
    Entrada: una matriz cuadrada A y un vector b.
    Salida: un vector x tal que Ax = b.
    """
    n = len(A)
    for k in range(n-1):
        if (A[k][k] != 0):
            M = matriz_de_eliminacion(A, k, 1)
            A = np.matmul(M, A)
            b = np.matmul(M, b)
        else:
            A, b = permutar(A, b, k)
    # FALTA:
    # 1. Permutación si el pivote es 0
    # 2. Decir que no hay solución si hay una columna llena de ceros
    x = sucesiva_hacia_atras(A, b)
    return x


# Método de eliminación de Gauss-Jordan
def gauss_jordan(A, b):
    """
    Entrada: una matriz cuadrada A y un vector b.
    Salida: un vector x tal que Ax = b.
    """
    n = len(A)
    for k in range(n):
        M = matriz_de_eliminacion(A, k, 0)
        A = np.matmul(M, A)
        b = np.matmul(M, b)
    x = [b[i] / A[i][i] for i in range(n)]
    return x


# Imprime la solución de un Sistema de Ecuaciones Lineales
def solucion_SEL(A, b, metodo):
    if (np.linalg.det(A) == 0):
        print("A es una matriz singular, el sistema no tiene solución.")
    else:
        if (metodo == 1):
            print("Método de Sucesión Sucesiva hacia atrás")
            inicio = time.time()
            x = sucesiva_hacia_atras(A, b)
            fin = time.time()
        elif (metodo == 2):
            print("Método de Sucesión Sucesiva hacia adelante")
            inicio = time.time()
            x = sucesiva_hacia_adelante(A, b)
            fin = time.time()
        elif (metodo == 3):
            print("Método de Eliminación de Gauss")
            inicio = time.time()
            x = gauss(A, b)
            fin = time.time()
        elif (metodo == 4):
            print("Método de Eliminación de Gauss-Jordan")
            inicio = time.time()
            x = gauss_jordan(A, b)
            fin = time.time()
        print("x = {0}\nTarda {1}\n".format(x, fin-inicio))


def main():

    # SUSTITUCIÓN SUCESIVA HACIA ATRÁS

    # Ejemplo 1
    A = [
         [6, -2, 2, 4],
         [0, -4, 2, 2],
         [0, 0, 2, -5],
         [0, 0, 0, -3],
    ]
    b = [12, 10, -9, -3]
    # solucion_SEL(A, b, 1)

    # Ejemplo 2
    A = [
         [4, 8, 3, 2, 9],
         [0, 1, 6, 3, 4],
         [0, 0, 9, 2, 5],
         [0, 0, 0, 7, 2],
         [0, 0, 0, 0, 1]
    ]
    b = [4, 7, 3, 8, 5]
    # solucion_SEL(A, b, 1)
    
    # Ejemplo 3
    A = [
         [7, 20, 13, 25, 89, 16],
         [0, 5, 56, 14, 77, 6],
         [0, 0, 17, 15, 10, 5],
         [0, 0, 0, 32, 8, 4],
         [0, 0, 0, 0, 9, 2],
         [0, 0, 0, 0, 0, 3]
    ]
    b = [8, 19, 50, 3, 24, 10]
    # solucion_SEL(A, b, 1)


    # SUSTITUCIÓN SUCESIVA HACIA ADELANTE
    # Ejemplo 1
    # Ejemplo 2
    # Ejemplo 3


    # GAUSS
    # Ejemplo 1
    A = [
        [5, 2, 7, 4],
        [2, 5, 1, 2],
        [8, 4, 6, -1],
        [-1, 2, -2, 6]
    ]
    b = [4, -2, 30, 21]
    # solucion_SEL(A, b, 3)

    # Ejemplo 2
    A = [
        [2, -1, 4, 1, -1],
        [-1, 3, -2, -1, 2],
        [5, 1, 3, -4, 1],
        [3, -2, -2, -3, 1],
        [-4, -1, -5, 3, -4]
    ]
    b = [7, 1, 33, 24, 49]
    # solucion_SEL(A, b, 3)

    # Ejemplo 3
    A = [
        [3, 10, 15, 45, 22, 97],
        [-13, 56, 70, -23, 60, 5],
        [8, 14, 13, -45, 41, 21],
        [73, -12, 54, 32, 61, 96],
        [12, -23, 25, 67, 19, 10],
        [-64, 32, 51, 73, -24, 34]
    ]
    b = [37, 71, 11, 42, 94, 96]
    # solucion_SEL(A, b, 3)
    print("x =",list(np.matmul(np.linalg.inv(A), b)))

    # GAUSS-JORDAN
    # Ejemplo 1
    # Ejemplo 2
    # Ejemplo 3


main()