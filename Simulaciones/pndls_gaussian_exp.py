
from procesos import *
from visualizacion import *
from Simulaciones.Input.inicializacion import *
from Simulaciones.Recursos.evolucion import *


if __name__ == '__main__':
    ####### INICIALIZACION #######
    eq = 'pndls_exp'
    bordes = 'periodic'
    fuente = 'gaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)
    sigma_forcing = (L / 3) / 2
    fuentes = fuente_pde(x_grid, Nx, Nt, source='gaussian', sigma=sigma_forcing)
    visualizacion(x_grid, t_grid, fuentes, tipo='colormap', guardar='no', path=simulation_data_path, file='ola',
                  nombre='plot')
    #plt.show()
    for i in range(1):
        for j in range(1):
            n = 4
            m = 1
            a = 1/3
            f = 14.8
            d = 20
            mu = 0.14

            ####### CONDICIONES INICIALES #########
            U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.001)
            V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.001)
            campos = campos_iniciales(Nt, Nx, [U_init, V_init])

            ####### DINAMICA #########
            campos_finales = RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, control=1, mu=mu, sigma=sigma_forcing, n=n, m=m, forcing_amp=a, forcing_freq=f, profundidad=d, forzamiento=fuentes)

            ####### DATOS #########
            campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)
            modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
            arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])
            sim_file = nombre_pndls_estandar(n=n, mu=mu, forcing_amp=a, forcing_freq=f, profundidad=d)
            if os.path.exists(simulation_data_path + sim_file) == True:
                print('Este archivo de ya existe, ¿desea eliminarlo y continuar? (y/n)')
                a = str(input())
                if a == 'y':
                    shutil.rmtree(simulation_data_path + sim_file)
                elif a == 'n':
                    sys.exit("Proceso terminado, cambie de carpeta")
            os.makedirs(simulation_data_path + sim_file)
            guardar_txt(simulation_data_path, sim_file, x_grid=x_grid, t_ligero=t_ligero, Phi_R=campo_ligeros[0], Phi_I=campo_ligeros[1])

            ####### VISUALIZACION #########
            visualizacion(x_grid, t_ligero, modulo, tipo='colormap', guardar='si', path=simulation_data_path, file=sim_file, nombre='plot')