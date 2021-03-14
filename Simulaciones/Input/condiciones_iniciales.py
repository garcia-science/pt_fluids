import numpy as np

def condiciones_iniciales_pde(type, x_grid, Nx, L, amplitude, params):
    if type == 'zero':
        U_init = amplitude * x_grid
    elif type == 'ones':
        U_init = amplitude * np.ones(Nx)
    elif type == 'sine':
        n = params[0]
        U_init = amplitude * np.sin(2 * np.pi * n * x_grid / L)
    elif type == 'cosine':
        n = params[0]
        U_init = amplitude * np.cos(2 * np.pi * n * x_grid / L)
    elif type == 'gaussian':
        if amplitude == 'norm':
            mu = params[0]
            sigma = params[1]
            U_init = (1 / sigma * np.sqrt(2 * np.pi))* np.exp(- (x_grid - mu) ** 2 / (2 * sigma ** 2))
        else:
            mu = params[0]
            sigma = params[1]
            U_init = amplitude * np.exp(- (x_grid - mu) ** 2/(2 * sigma ** 2))
    elif type == 'soliton':
        c = params[0]
        U_init = 3 * c * (1 / np.cosh(np.sqrt(c / 4) * x_grid)) ** 2
    elif type == '2soliton':
        c_1 = params[0]
        c_2 = c_1 / 2
        centro_1 = params[1]
        centro_2 = centro_1 + params[2]
        U_init_1 = 3 * c_1 * (1 / np.cosh(np.sqrt(c_1 / 4) * (x_grid - centro_1))) ** 2
        U_init_2 = + 3 * c_2 * (1 / np.cosh(np.sqrt(c_2 / 4) * (x_grid - centro_2))) ** 2
        U_init = U_init_1 + U_init_2
    return U_init
