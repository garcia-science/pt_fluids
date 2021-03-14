from Experimento.Visualizacion.tres_dimensiones import *
from Source.RK4 import *
from matplotlib import pyplot as plt
from Source.preparacion import *
from Input.condiciones_iniciales import *
from Input.fuentes import *
from Input.parametros import *
from matplotlib.animation import FuncAnimation


if __name__ == '__main__':

    eq = 'KdV'
    bordes = 'periodic'
    fuente = 'gaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)

    ####### PARAMETROS Y FUENTES #########
    fuentes = fuente_pde(x_grid, t_grid, Nx, Nt, 'none')
    #fig_1 = color_map(x_grid, t_grid, fuentes)
    #plt.show()
    parametros = parametros_PDE(eq, fuentes)

    ####### CONDICIONES INICIALES #########
    U_init = condiciones_iniciales_pde('2soliton', x_grid, Nx, L, 1, [1, -15, 10])
    V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.001, [0, 0])
    condiciones_iniciales = [U_init, V_init]
    campos = campos_iniciales(Nt, Nx, condiciones_iniciales)

    ####### DINAMICA #########
    campos_finales = RK4_PDE(campos, bordes, parametros, dx, dt, Nx, Nt, eq, x_grid, t_grid)
    #A = campos_finales[0]**2 + campos_finales[1]**2
    #fig_1 = color_map(x_grid, t_grid, campos_finales[0])
    #plt.show()

    x = x_grid

    fig = plt.figure()

    lines = plt.plot([])
    line = lines[0]
    #other setuo
    plt.xlim(-L/2, L/2)
    plt.ylim(-1.1, 10)
    def animate(frame):
        y = campos_finales[0][100 * frame, :]
        line.set_data((x, y))
        #update plot

    anim = FuncAnimation(fig, animate, frames=int(Nt/100), interval=20)
    anim.save('basic_animation.gif', fps=100)
    #video = anim.to_html5_video()
    #html = display.HTML(video)
    #display.display(html)
    plt.show()












