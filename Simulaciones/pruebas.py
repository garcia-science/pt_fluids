from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
import numpy as np
from matplotlib import pyplot as plt
from Source.preparacion import *

if __name__ == '__main__':
    #Parametros de grilla
    dx = 0.01
    dt = 0.01
    c = 1
    C = c*(dt/dx)
    x_min = -10
    x_max = 10
    L = x_max - x_min
    T = 20
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    print(Nx)
    print(Nt)
    u_min = u_max = 0


    #I_x = np.sin(2 * np.pi * x_grid / L)
    I_x = 100*np.exp(- x_grid ** 2)
    U = I_x.tolist()
    U_final = [U]
    A = [0]*int(Nx)
    for j in range(int(Nx) - 1):
        if j == 1:
            A[j] = U[j] - 0.5 * C ** 2 * (U[j + 1] - 2 * U[j] + u_min)
        elif j == int(Nx) - 2:
            A[j] = U[j] - 0.5 * C ** 2 * (u_max - 2 * U[j] + U[j - 1])
        elif j == 0:
            A[j] = u_min
        elif j == int(Nx) - 1:
            A[j] = u_max
        else:
            A[j] = U[j] - 0.5 * C ** 2 * (U[j + 1] - 2 * U[j] + U[j - 1])
    U_final.append(A)
    print(U_final)
    for i in range(1, int(Nt) - 1):
        j = 0
        while j <= int(Nx) - 1:
            #print(str(i)+", "+str(j))
            if j == 1:
                A[j] = -1 * U_final[i - 1][j] + 2 * U_final[i][j] + C ** 2 * (
                            U_final[i][j + 1] - 2 * U_final[i][j] + u_min)
            elif j == int(Nx) - 2:
                A[j] = -1 * U_final[i - 1][j] + 2 * U_final[i][j] + C ** 2 * (
                            u_max - 2 * U_final[i][j] + U_final[i][j - 1])
            elif j == 0:
                A[j] = u_min
            elif j == int(Nx) - 1:
                A[j] = u_max
            else:
                A[j] = -1 * U_final[i - 1][j] + 2 * U_final[i][j] + C ** 2 * (
                        U_final[i][j + 1] - 2 * U_final[i][j] + U_final[i][j - 1])
            j += 1
        U_final.append(A)
        #print(U_final[i][:])
    U_f = np.array(U_final)
    print(np.shape(U_f))
    fig = color_map(x_grid, t_grid, U_f)
    plt.show()







    #PSI_inicial = campo_inicial(Nx_pasos, Nt_pasos, psi_0_2)

