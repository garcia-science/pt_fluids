import numpy as np
from procesos import *
from visualizacion import *
from Simulaciones.Input.inicializacion import *
from Simulaciones.Recursos.evolucion import *

if __name__ == '__main__':
    eq = 'pndls'
    bordes = 'periodic'
    fuente = 'gaussian'
    ################---ADIMENSIONAL---#####################
    #L = 480
    #sigma_forcing = L/3
    #alpha = 1
    #beta = 1
    #gamma = 0.225
    #mu = 0.15
    #nu = 0.15
    ################---DIMENSIONAL---#####################
    unidad = 0.1
    g = 9790 * unidad
    L = 480 * unidad
    l_x = L / 3
    l_y = 16 * unidad
    d = 20 * unidad
    n = 5
    m = 1
    a_ang = 7
    f = 15.00

    a_mm = (28 + a_ang / 12) * unidad
    sigma_forcing = l_x / (2 * 2.3482)
    k_x = np.pi * n / l_x
    k_y = np.pi * m / l_y
    k = np.sqrt(k_x ** 2 + k_y ** 2)
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau) / (2 * np.pi)
    w = f / 2
    GAMMA = 4 * w ** 2 * a_mm

    alpha = (1 / (4 * k ** 2)) * (1 + k * d * ((1 - tau ** 2) / tau)) #término difusivo
    beta = (k ** 2 / 64) * (6 * tau ** 2 - 5 + 16 * tau ** (-2) - 9 * tau ** (-4)) #término no lineal
    gamma = GAMMA / (4 * g) #amplitud adimensional
    nu = 0.5 * ((w / w_1) ** 2 - 1)
    mu = 0.0125

    ################---DEFINICIÓN DE LA MALLA---#####################
    [xmin, xmax, dx] = [-L/2, L/2, 0.2]
    [tmin, tmax, dt] = [0, 6000, 0.01]
    x_grid = np.arange(xmin, xmax + dx, dx)
    t_grid = np.arange(tmin, tmax + dt, dt)
    T = tmax
    Nx = x_grid.shape[0]
    Nt = t_grid.shape[0]
    fuentes = fuente_pde(x_grid, Nx, Nt, source='gaussian', sigma=sigma_forcing, dist=0)
    fuente_ligera, t_ligero = campos_ligeros([fuentes, fuentes], 100, Nt, Nx, T)
    visualizacion(x_grid, t_ligero, fuente_ligera[0], tipo='colormap', guardar='no', path=simulation_data_path,
                  file="a", nombre='plot')
    plt.show()
    ####### CONDICIONES INICIALES #########
    U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.01)
    V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.1)
    campos = campos_iniciales(Nt, Nx, [U_init, V_init])

    ####### DINAMICA #########
    campos_finales = RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, control=0.01, alpha=alpha, beta=beta, gamma=gamma,
                             mu=mu, nu=nu, forzamiento=fuentes)

    ####### DATOS #########
    campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)
    modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
    arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])
    sim_file = nombre_pndls_estandar(gamma=gamma, mu=mu, nu=nu)
    if os.path.exists(simulation_data_path + sim_file):
        shutil.rmtree(simulation_data_path + sim_file)
    os.makedirs(simulation_data_path + sim_file)
    guardar_txt(simulation_data_path, sim_file, X=x_grid, T=t_ligero, real=campo_ligeros[0],
                img=campo_ligeros[1], mod=modulo, arg=arg)
    ####### VISUALIZACION #########
    visualizacion(x_grid, t_ligero, campo_ligeros[0], tipo='colormap', guardar='si', path=simulation_data_path,
                  file=sim_file, nombre='plot')
    plt.show()