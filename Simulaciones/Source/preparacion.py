import numpy as np

def grilla(x_min, x_max, T, dx, dt):
    Nx_pasos = (x_max - x_min) / dx
    Nt_pasos = T / dt
    print('tenemos ' + str(Nx_pasos) + ' pasos espaciales y ' + str(Nt_pasos) + ' pasos temporales')
    x_grid = np.linspace(x_min, x_max, int(Nx_pasos))
    t_grid = np.linspace(0, T, int(Nt_pasos))
    return Nx_pasos, Nt_pasos, x_grid, t_grid
def campo_inicial(Nx_pasos, Nt_pasos, phi_inicial):
    PHI = np.zeros((int(Nt_pasos) - 1, int(Nx_pasos)))
    PHI = np.append([phi_inicial], PHI, axis=0)
    return PHI