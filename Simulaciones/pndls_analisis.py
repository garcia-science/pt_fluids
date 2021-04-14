import os
from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
from matplotlib import pyplot as plt
from Source.preparacion import *
from Input.condiciones_iniciales import *
from Input.fuentes import *

if __name__ == '__main__':
    gamma = str(0.2)
    mu = str(0.15)
    nu = str(-0.2)
    file = nombre_pndls(gamma, mu, nu)
    [x_grid, t_grid, U_r, U_i] = cargar_txt(4, sim_path, file)
    print('Datos extraidos')
    modulo = np.sqrt(U_r ** 2 + U_i ** 2)
    arg = np.arctan2(U_i, U_r)
    A = modulo
    #B = A.real
    t_grid = np.linspace(0, 400, int(400 / (100 * 0.001)))
    print('Módulo y fase calculados')
    color_map_for(x_grid, t_grid, arg, guardar='no', path=sim_path, file=file, nombre='/mod',
                  titulo='$\gamma = $' + gamma +'$\ \mu = $' + mu + '$\ \\nu = $' + nu, xlabel=r'$x$ (Espacio)', ylabel=r'$t$ (tiempo)',
                  z_name=r'       $|\psi(x,t)|$')
    plt.show()