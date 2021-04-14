import os
from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
from matplotlib import pyplot as plt
from Source.preparacion import *
from Input.condiciones_iniciales import *
from Input.fuentes import *
from Experimento.Source.procesos import *

if __name__ == '__main__':
    ####### INICIALIZACION #######
    eq = 'pndls'
    bordes = 'periodic'
    fuente = 'gaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)

    ####### PARAMETROS Y FUENTES #########
    fuentes = fuente_pde(x_grid, t_grid, Nx, Nt,  'bigaussian')

    for i in range(1):
        for j in range(1):
            parametros = [fuentes, 1, 1, 1, 1 + j * 0.05, 0.1, 1 + 0.05 * i]
            #delta = - parametros[6] + np.sqrt(parametros[4] ** 2 - parametros[5] ** 2)
            gamma = str(round(parametros[4], 3))
            mu = str(round(parametros[5], 3))
            nu = str(round(parametros[6], 3))
            ####### CONDICIONES INICIALES #########
            U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.01, [0, 3])
            V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.01, [0, 3])
            campos = campos_iniciales(Nt, Nx, [U_init, V_init])

            ####### DINAMICA #########
            campos_finales = RK4_PDE(campos, bordes, parametros, dx, dt, Nx, Nt, eq, x_grid, t_grid)

            ####### DATOS #########
            campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)
            modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
            arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])
            file = nombre_pndls_bigaussian(gamma, mu, nu, sigma, dist)
            os.mkdir(sim_path + file)
            guardar_txt([x_grid, t_ligero, campo_ligeros[0], campo_ligeros[1]], sim_path, file)

            ####### VISUALIZACION #########
            color_map_for(x_grid, t_ligero, modulo, guardar='si', path=sim_path, file=file, nombre='/mod',
                          titulo='$\gamma = $' + gamma + '$\ \mu = $' + mu + '$\ \\nu = $' + nu, xlabel=r'$x$ (Espacio)',
                          ylabel=r'$t$ (tiempo)',
                          z_name=r'          $|\psi(x,t)|$')
            plt.close()
