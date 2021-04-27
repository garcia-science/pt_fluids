
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
    sigma_forcing = (L / 3) / 2
    fuentes = fuente_pde(x_grid, t_grid, Nx, Nt, sigma=sigma_forcing, dist=0, source='gaussian')

    for i in range(1):
        for j in range(1):
            n = 4
            a = 28.32
            w = 0.62
            d = 2
            mu = 0.03

            ####### CONDICIONES INICIALES #########
            U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.001)
            V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.001)
            campos = campos_iniciales(Nt, Nx, [U_init, V_init])

            ####### DINAMICA #########
            campos_finales = RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, control=1, mu=mu, sigma=sigma_forcing, n=n, forcing_amp=a, forcing_freq=w, profundidad=d, forzamiento=fuentes)

            ####### DATOS #########
            campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)
            modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
            arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])
            sim_file = nombre_pndls_estandar(n=n, mu=mu, forcing_amp=a, forcing_freq=w, profundidad=d)
            if os.path.exists(simulation_data_path + sim_file) == True:
                print('Este archivo de CANNY ya existe, ¿desea eliminarlo y continuar? (y/n)')
                a = str(input())
                if a == 'y':
                    shutil.rmtree(simulation_data_path + sim_file)
                elif a == 'n':
                    sys.exit("Proceso terminado, cambie de carpeta")
            os.makedirs(simulation_data_path + sim_file)
            guardar_txt(simulation_data_path, sim_file, x_grid=x_grid, t_ligero=t_ligero, Phi_R=campo_ligeros[0], Phi_I=campo_ligeros[1])

            ####### VISUALIZACION #########
            visualizacion(x_grid, t_ligero, modulo, tipo='colormap', guardar='si', path=simulation_data_path, file=sim_file, nombre='plot')