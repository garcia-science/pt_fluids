import numpy as np


def grilla(x_min, x_max, t, dx, dt):
    Nx_pasos = (x_max - x_min) / dx
    Nt_pasos = t / dt
    print('tenemos ' + str(Nx_pasos) + ' pasos espaciales y ' + str(Nt_pasos) + ' pasos temporales')
    x_grid = np.linspace(x_min, x_max, int(Nx_pasos))
    t_grid = np.linspace(0, t, int(Nt_pasos))
    return Nx_pasos, Nt_pasos, x_grid, t_grid


def campos_iniciales(nt, nx, condiciones_iniciales):
    N_fields = len(condiciones_iniciales)
    campos = [0] * N_fields
    for i in range(N_fields):
        U = np.zeros((nt, nx))
        U[0, :] = condiciones_iniciales[i]
        campos[i] = U
    return campos
