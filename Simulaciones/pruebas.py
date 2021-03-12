from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
import numpy as np
from matplotlib import pyplot as plt
from Source.preparacion import *

if __name__ == '__main__':
    #Parametros de grilla
    dx = 1
    dt = 0.001
    x_min = -200
    x_max = 200
    L = x_max - x_min
    T = 4
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    Nx = int(Nx)
    Nt = int(Nt)

    eq = 'pndls'
    bordes = 'fixed'
    sigma = 4
    frec = 0.75
    perfil_forzamiento = np.exp(- x_grid ** 2 / (2 * sigma ** 2))
    forzamiento = [perfil_forzamiento] * Nt
    forzamiento = np.array(forzamiento)
    forzamiento = np.transpose(forzamiento)
    for j in range(Nx):
        forzamiento[j, :] = forzamiento[j, :] * np.sin(frec * t_grid)
    forzamiento = np.transpose(forzamiento)
    fig_1 = color_map(x_grid, t_grid, forzamiento)
    plt.show()
    parametros = [forzamiento, 1, 0.1, 0.01, 0.1]  # que se generen solo con la ecuacion

    '''
    CONDICIONES INICIALES (hacer set de condiciones usuales y modular)
    '''
    #U_init = 0.01 * np.sin(2 * np.pi * x_grid / L)
    #U_init = 0.0001*np.random.rand(1, Nx)
    #U_init = 0 * x_grid #posicion cero
    U_init = 0.001 * np.exp(- x_grid ** 2/30) #posicion gaussiana
    #U_init = (1 / np.cosh(x_grid)) ** 2

    V_init = 0 * x_grid#velocidad cero
    #V_init = np.exp(- x_grid ** 2 / (2 * 10))  # posicion gaussiana
    #sigma = 1
    #V_init = (1 / np.sqrt(2) * np.pi * sigma) * np.exp(- x_grid ** 2/sigma) #velocidad gaussiana

    condiciones_iniciales = [U_init, V_init]
    campos = campos_iniciales(Nt, Nx, condiciones_iniciales)

    '''
    DINAMICA
    '''

    campos_finales = RK4_PDE(campos, bordes, parametros, dx, dt, Nx, Nt, eq, x_grid, t_grid)
    A = campos_finales[0]**2 + campos_finales[1]**2
    fig_1 = color_map(x_grid, t_grid, A)
    plt.show()












