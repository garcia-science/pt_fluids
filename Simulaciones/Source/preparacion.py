import numpy as np


def iniciar_PDE(eq):
    if eq == 'pndls':
        dx = 0.5
        dt = 0.001
        x_min = -30
        x_max = 30
        l = x_max - x_min
        t = 1000
    elif eq == 'wave':
        dx = 0.5
        dt = 0.0005
        x_min = -50
        x_max = 50
        l = x_max - x_min
        t = 10
    elif eq == 'fisher':
        dx = 0.5
        dt = 0.005
        x_min = -50
        x_max = 50
        l = x_max - x_min
        t = 30
    elif eq == 'KdV':
        dx = 0.25
        dt = 0.0005
        x_min = -30
        x_max = 30
        l = x_max - x_min
        t = 30
    return dx, dt, x_min, x_max, l, t


def grilla(x_min, x_max, t, dx, dt):
    Nx_pasos = (x_max - x_min) / dx
    Nt_pasos = t / dt
    print('tenemos ' + str(Nx_pasos) + ' pasos espaciales y ' + str(Nt_pasos) + ' pasos temporales')
    x_grid = np.linspace(x_min, x_max, int(Nx_pasos))
    t_grid = np.linspace(0, t, int(Nt_pasos))
    Nx_pasos = int(Nx_pasos)
    Nt_pasos = int(Nt_pasos)
    return Nx_pasos, Nt_pasos, x_grid, t_grid


def campos_iniciales(nt, nx, condiciones_iniciales):
    N_fields = len(condiciones_iniciales)
    campos = [0] * N_fields
    for i in range(N_fields):
        U = np.zeros((nt, nx))
        U[0, :] = condiciones_iniciales[i]
        campos[i] = U
    return campos
