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
    fuente = 'bigaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)

    ####### PARAMETROS Y FUENTES #########
    for i in range(1):
        sigma = 3
        sigma_st = str(sigma)
        dist = 40 + i
        dist_st = str(dist)
        fuentes = fuente_pde(x_grid, t_grid, Nx, Nt, sigma, dist, fuente)
        parametros = [fuentes, 1, 1, 1, 0.2, 0.1, 0.2]
        gamma = str(round(parametros[4], 3))
        mu = str(round(parametros[5], 3))
        nu = str(round(parametros[6], 3))
        file = nombre_pndls_bigaussian(gamma, mu, nu, sigma_st, dist_st)
        os.makedirs(sim_path + file)
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
        guardar_txt([x_grid, t_ligero, campo_ligeros[0], campo_ligeros[1]], sim_path, file)

        ####### VISUALIZACION #########
        color_map_for(x_grid, t_ligero, modulo, guardar='si', path=sim_path, file=file, nombre='/mod',
                      titulo='$\sigma =$ ' + sigma_st + '$d =$ ' + dist_st, xlabel=r'$x$ (Espacio)',
                      ylabel=r'$t$ (tiempo)',
                      z_name=r'          $|\psi(x,t)|$')
        plt.close()
