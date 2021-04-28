
from procesos import *
from visualizacion import *
from Simulaciones.Input.inicializacion import *
from Simulaciones.Recursos.evolucion import *

if __name__ == '__main__':
    eq = 'pndls'
    bordes = 'periodic'
    fuente = 'bigaussian'
    dx, dt, x_min, x_max, L, T = iniciar_PDE(eq)
    Nx, Nt, x_grid, t_grid = grilla(x_min, x_max, T, dx, dt)

    for i in range(1):
        for j in range(1):
            gamma = 0.28
            mu = 0.1
            nu = 0.32

            sigma = 3 + i
            distancia = 20 + j
            fase = np.pi
            fuentes = fuente_pde(x_grid, Nx, Nt, source=fuente, sigma=sigma, distancia=distancia, fase=fase)

            U_init = condiciones_iniciales_pde('ones', x_grid, Nx, L, 0.01)
            V_init = condiciones_iniciales_pde('zero', x_grid, Nx, L, 0.01)
            campos = campos_iniciales(Nt, Nx, [U_init, V_init])

            campos_finales = RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, control=1, gamma=gamma, mu=mu, nu=nu,
                                     forzamiento=fuentes)
            campo_ligeros, t_ligero = campos_ligeros(campos_finales, 100, Nt, Nx, T)

            modulo = np.sqrt(campo_ligeros[0] ** 2 + campo_ligeros[1] ** 2)
            arg = np.arctan2(campo_ligeros[0], campo_ligeros[1])

            sim_file = nombre_pndls_bigaussian(gamma, mu, nu, sigma, distancia, fase)
            #if os.path.exists(simulation_data_path + '\\mu=' + str(mu) + sim_file) == True:
            #    print('Este directorio ya existe, ¿deseas eliminarlo y continuar? (Y/N)')
            #    a = str(input())
            #    if a == 'Y' or 'y':
            #        shutil.rmtree(simulation_data_path + sim_file)
            #    elif a == 'N' or 'n':
            #        sys.exit("Proceso terminado, cambie de carpeta")
            os.makedirs(simulation_data_path + sim_file)
            guardar_txt(simulation_data_path, sim_file, x_grid=x_grid, t_grid=t_ligero,
                        PHI_real=campo_ligeros[0], PHI_img=campo_ligeros[1], PHI_modulo=modulo, PHI_arg=arg)

            visualizacion(x_grid, t_ligero, modulo, tipo='colormap', guardar='si', path=simulation_data_path, file=sim_file, nombre='plot')
