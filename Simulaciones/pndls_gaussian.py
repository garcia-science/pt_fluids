
from procesos import *
from visualizacion import *
from Simulaciones.Input.inicializacion import *
from Simulaciones.Recursos.evolucion import *


if __name__ == '__main__':
    ####### INICIALIZACION #######
    eq = 'pndls'
    bordes = 'periodic'
    fuente = 'gaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    largo_celda = 480
    sigma_forcing = L/3
    fuentes = fuente_pde(x_grid, t_grid, Nx, Nt, sigma=sigma_forcing, dist=0, source='gaussian')

    for i in range(1):
        for j in range(1):
            alpha = 1
            beta = 1
            gamma = 0.225
            mu = 0.15 + 0.1 * i
            nu = 0.15 + 0.1 * j

            ####### CONDICIONES INICIALES #########
            U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.01)
            V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.01)
            campos = campos_iniciales(Nt, Nx, [U_init, V_init])

            ####### DINAMICA #########
            campos_finales = RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, control=1, alpha=alpha, beta=beta, gamma=gamma,
                                     mu=mu, nu=nu, forzamiento=fuentes)

            ####### DATOS #########
            campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)
            modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
            arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])
            sim_file = nombre_pndls_estandar(gamma=gamma, mu=mu, nu=nu)
            os.makedirs(simulation_data_path + sim_file)
            guardar_txt(simulation_data_path, sim_file, X=x_grid, T=t_ligero, real=campo_ligeros[0],
                        img=campo_ligeros[1], mod=modulo, arg=arg)
            ####### VISUALIZACION #########
            visualizacion(x_grid, t_ligero, modulo, tipo='colormap', guardar='si', path=simulation_data_path, file=sim_file, nombre='plot')
