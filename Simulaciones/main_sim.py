from Input.parametros import parametros
from Input.variables import *
from Source.ODE import *
from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
import numpy as np
from matplotlib import pyplot as plt
from Source.preparacion import *

if __name__ == '__main__':
############1D EQUATION ODE, FIRST ORDER############
    #eqn, params, y_0 = parametros()
    #h, T, n, x = variables(y_0)
    #f = funcion(eqn, x, params)
    #y_t = RK4_ODE(f, x, h, T, eqn, params)
    #plt.plot(T, y_t)
    #plt.show()

############1D EQUATION PDE, FIRST ORDER############
    #Parametros de grilla
    dx = 0.1
    dt = 0.01
    x_min = -10
    x_max = 10
    T = 100

    #Parametros
    eq = 'KdV'
    c_1 = 0.01
    c_2 = 0.001
    parametros = [c_1, c_2]


    #Preparacion de grilla, condicion inicial y condiciones de borde
    Nx_pasos, Nt_pasos, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    BC = 'periodicas'
    phi_inicial = (1 / np.cosh(x_grid))**2
    PHI_inicial = campo_inicial(Nx_pasos, Nt_pasos, phi_inicial)


    #EVOLUCIÓN TEMPORAL
    PHI = RK4_PDE(PHI_inicial, dx, dt, Nt_pasos, Nx_pasos, eq, parametros, BC)

    #PLOTEO
    #XX, TT = np.meshgrid(x_grid, t_grid)
    #fig = plot_ZXT(x_grid, t_grid, PHI)
    fig = color_map(x_grid, t_grid, PHI_inicial)
    #plt.plot(T, y_t)
    plt.show()
