from Input.variables import *
from Source.ODE import *
from Experimento.Visualizacion.tres_dimensiones import *
from Experimento.Input.enrutador import *
from Source.RK4 import *
from matplotlib import pyplot as plt
from Source.preparacion import *
from Input.condiciones_iniciales import *
from Input.fuentes import *
from Input.parametros import *
from matplotlib.animation import FuncAnimation

if __name__ == '__main__':
############1D EQUATION ODE, FIRST ORDER############
    #eqn, params, y_0 = parametros()
    #h, T, n, x = variables(y_0)
    #f = funcion(eqn, x, params)
    #y_t = RK4_ODE(f, x, h, T, eqn, params)
    #plt.plot(T, y_t)
    #plt.show()

############1D EQUATION PDE############
    ####### INICIALIZACION #######
    eq = 'pndls'
    bordes = 'periodic'
    fuente = 'gaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)

    ####### PARAMETROS Y FUENTES #########
    fuentes = fuente_pde(x_grid, t_grid, Nx, Nt, 'gaussian')
    fig_1 = color_map(x_grid, t_grid, fuentes, guardar = 'no', nombre = 'no', titulo = 'Perfil de Forzamiento', xlabel=r'$x$ (Espacio)', ylabel=r'$t$ (tiempo)')
    plt.show()
    parametros = parametros_PDE(eq, fuentes)

    ####### CONDICIONES INICIALES #########
    U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.001, [1, -15, 10])
    V_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.001, [0, 0])
    campos = campos_iniciales(Nt, Nx, [U_init, V_init])

    ####### DINAMICA #########
    campos_finales = RK4_PDE(campos, bordes, parametros, dx, dt, Nx, Nt, eq, x_grid, t_grid)
    A = campos_finales[0]**2 + campos_finales[1]**2

    ####### VISUALIZACION #########
    fig_1 = color_map(x_grid, t_grid, A, guardar = 'no', nombre = 'no', titulo = 'Diagrama espacio-temporal ', xlabel=r'$x$ (Espacio)', ylabel=r'$t$ (tiempo)')
    #anim = animacion(campos_finales[0], x_grid, Nt, L, [-2, 2], guardar='si', nombre='/anim_01', titulo='Animación PDE', xlabel=r'$x$ (Espacio)', ylabel=r'$t$ (tiempo)')
    plt.show()