from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    #Parametros de grilla
    dx = 0.1
    dt = 0.01
    x_min = -10
    x_max = 10
    T = 100
    Nx_pasos = (x_max - x_min)/dx
    Nt_pasos = T/dt
    print('tenemos '+str(Nx_pasos)+' pasos espaciales y '+str(Nt_pasos)+' pasos temporales')

    #Parametros ecuacion de fisher
    eq = 'KdV'
    D = 0.01
    r = 0.001
    parametros = [D, r]


    #Preparacion de grilla, condicion inicial y condiciones de borde
    BC = 'periodicas'
    x_grid = np.linspace(x_min, x_max, int(Nx_pasos))
    t_grid = np.linspace(0, T, int(Nt_pasos))
    phi_inicial = 2*np.exp(-x_grid**2)

    PHI = np.zeros((int(Nx_pasos), int(Nt_pasos)-1))
    PHI = np.matrix.transpose(PHI)
    PHI = np.append([phi_inicial], PHI, axis=0)
    print(PHI.shape)
    print(phi_inicial.shape)
    print(t_grid.shape)

    PHI = RK4_PDE(PHI, dx, dt, Nt_pasos, Nx_pasos, eq, parametros, BC)

    #PLOTEO
    XX, TT = np.meshgrid(x_grid, t_grid)
    #fig = plot_ZXT(XX, TT, PHI)
    fig = color_map(x_grid, t_grid, PHI)
    plt.show()





