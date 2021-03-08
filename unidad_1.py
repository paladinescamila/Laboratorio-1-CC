# UNIDAD 1: SISTEMAS DE PUNTO FLOTANTE

import matplotlib.pyplot as plt
import time


# Conversor de entero a binario
def calcular_di(i, t):
    """
    Entrada: dos enteros i y j.
    Salida: cadena de la representación binaria de i, con una longitud t.
    """
    ceros = ""
    binario = bin(i)[2:]
    for j in range(t - 1 - len(binario)):
        ceros += "0"

    return ceros + binario


# Generador de números del Sistema de Punto Flotante (con beta = 2)
def punto_flotante(t, L, U):
    """
    Entrada: tres enteros t, L y U.
    Salida: lista de los números del sistema de punto flotante, en base 2, de 
            precisión t y [L, U] como rango del exponente.
    """
    numeros = [0]    
    for e in range(L, U + 1):
        for i in range(t + 1):
            x = 1
            di = calcular_di(i, t)
            for j in range(t - 1):
                x += int(di[j]) / (2**(j+1))            
            numeros += [x*(2**e), x*(2**e)*(-1)]
    
    numeros = sorted(list(set(numeros)))
    N = (2**t) * (U-L+1) + 1
    UFL = 2**L
    OFL = 2**(U+1) * (1-2**(-t))
    
    return numeros, N, UFL, OFL


# Imprime y muestra la gráfica para un Sistema de Punto Flotante
def dibujar_SPF(t, L, U):
    inicio = time.time()
    numeros, N, UFL, OFL = punto_flotante(t, L, U)
    fin = time.time()

    print("Sistema: {0}".format(numeros))
    print("N: {0}\nUFL: {1}\nOFL: {2}".format(N, UFL, OFL))
    print("Tarda {0}".format(fin-inicio))

    for i in numeros:
        plt.plot(i, 0, marker="o", color="blue", markersize = 10)
    plt.grid()
    plt.show()
    print()


# EJEMPLOS DE PRUEBA (También se encuentran en el informe)
def main():
	
    # Ejemplo 1
    dibujar_SPF(3,-1,3)

    # Ejemplo 2
    dibujar_SPF(2,-2,2)

    # Ejemplo 3
    dibujar_SPF(3,-2,1)


main()