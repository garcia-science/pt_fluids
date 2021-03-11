from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
import numpy as np
from matplotlib import pyplot as plt
from Source.preparacion import *

if __name__ == '__main__':
    #Parametros de grilla
    dx = 0.5
    dt = 0.01
    c = 1
    x_min = -50
    x_max = 50
    L = x_max - x_min
    T = 200
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    Nx = int(Nx)
    Nt = int(Nt)
    U = np.zeros((Nx, Nt))
    V = np.zeros((Nx, Nt))
    U_init = np.sin(2 * 2 * np.pi * x_grid / L) #posicion sinosidal
    U_init = 0 * x_grid #posicion cero
    #U_init = np.exp(- x_grid ** 2) #posicion gaussiana

    #V_init = 0 * x_grid #velocidad cero
    V_init = np.exp(- x_grid ** 2) #velocidad gaussiana

    U[:, 0] = U_init
    V[:, 0] = V_init
    U_left = U_right = 0
    V_left = V_right = 0
    def funciones(i, j, U, V, Nx, U_left, U_right, V_left, V_right, c):
        if i == 0:
            F = V_left
            G = 0
        elif i == 1:
            F = V[i, j]
            G = (c ** 2) * (U[i + 1, j] - 2 * U[i, j] + U_left)
        elif i == Nx - 1:
            F = V_right
            G = 0
        elif i == Nx - 2:
            F = V[i, j]
            G = (c ** 2) * (U_right - 2 * U[i, j] + U[i - 1, j])
        else:
            F = V[i, j]
            G = (c ** 2) * (U[i + 1, j] - 2 * U[i, j] + U[i - 1, j])
        return F, G

    for j in range(0, Nt - 1):
        print('dt = ' + str(j))
        for i in range(0, Nx - 1):
            F, G = funciones(i, j, U, V, Nx, U_left, U_right, V_left, V_right, c)
            k1 = F
            m1 = G
            k2 = F + 0.5 * dt * k1
            m2 = G + 0.5 * dt * k1
            k3 = F + 0.5 * dt * k2
            m3 = G + 0.5 * dt * k2
            k4 = F + dt * k3
            m4 = G + dt * k3

            V[i, j + 1] = V[i, j] + dt * (m1 + 2 * m2 + 2 * m3 + m4) / 6
            U[i, j + 1] = U[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    fig_1 = color_map(x_grid, t_grid, np.transpose(U))
    #fig_2 = color_map(x_grid, t_grid, np.transpose(V))
    plt.show()












