import numpy as np


def iniciar_PDE(eq):
    if eq == 'pndls':
        dx = 0.5
        dt = 0.001
        x_min = -50
        x_max = 50
        l = x_max - x_min
        t = 600
    elif eq == 'pndls_exp':
        dx = 0.5
        dt = 0.001
        x_min = -240
        x_max = 240
        l = x_max - x_min
        t = 5
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
    print('Hay ' + str(Nx_pasos) + ' pasos espaciales y ' + str(Nt_pasos) + ' pasos temporales')
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


def fuente_pde(x_grid, Nx, Nt, source, **kwargs):
    if source == 'gaussian':
        sigma = kwargs['sigma']
        fuente_init = np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        fuente = [fuente_init] * Nt
        fuente = np.array(fuente)
    elif source == 'bigaussian':
        sigma = kwargs['sigma']
        distancia = kwargs['distancia']
        phase = kwargs['fase']
        fuente_init_1 = np.exp(- (x_grid + distancia / 2) ** 2 / (2 * sigma ** 2))
        fuente_init_2 = np.cos(phase) * np.exp(- (x_grid - distancia / 2) ** 2 / (2 * sigma ** 2))
        fuente_1 = np.array([fuente_init_1] * Nt)
        fuente_2 = np.array([fuente_init_2] * Nt)
        fuente = fuente_1 + fuente_2
    elif source == 'none':
        fuente = np.array([[0] * Nx] * Nt)
    elif source == 'pistons':
        print('aun no me programas')
    return fuente


def condiciones_iniciales_pde(tipo, x_grid, Nx, L, amplitude, **kwargs):
    if tipo == 'zero':
        U_init = amplitude * np.zeros(Nx)
    elif tipo == 'ones':
        U_init = amplitude * np.ones(Nx)
    elif tipo == 'sine':
        n = kwargs['n']
        U_init = amplitude * np.sin(2 * np.pi * n * x_grid / L)
    elif tipo == 'cosine':
        n = kwargs['n']
        U_init = amplitude * np.cos(2 * np.pi * n * x_grid / L)
    elif tipo == 'gaussian':
        if amplitude == 'norm':
            mu = kwargs['mu']
            sigma = kwargs['sigma']
            U_init = (1 / sigma * np.sqrt(2 * np.pi))* np.exp(- (x_grid - mu) ** 2 / (2 * sigma ** 2))
        else:
            mu = kwargs['mu']
            sigma = kwargs['sigma']
            U_init = amplitude * np.exp(- (x_grid - mu) ** 2/(2 * sigma ** 2))
    elif tipo == 'soliton':
        c = kwargs['c']
        U_init = 3 * c * (1 / np.cosh(np.sqrt(c / 4) * x_grid)) ** 2
    elif tipo == '2soliton':
        c_1 = kwargs['velocidad']
        c_2 = c_1 / 2
        centro_1 = kwargs['posicion_1']
        centro_2 = kwargs['posicion_2']
        U_init_1 = 3 * c_1 * (1 / np.cosh(np.sqrt(c_1 / 4) * (x_grid - centro_1))) ** 2
        U_init_2 = + 3 * c_2 * (1 / np.cosh(np.sqrt(c_2 / 4) * (x_grid - centro_2))) ** 2
        U_init = U_init_1 + U_init_2
    elif tipo == 'soliton_pndls':
        nu = kwargs['nu']
        gamma = kwargs['gamma']
        mu = kwargs['mu']
        delta = - nu + np.sqrt(gamma ** 2 - mu ** 2)
        U_init = np.sqrt(2 * delta) * (1 / np.cosh(np.sqrt(delta) * x_grid))
    return U_init

