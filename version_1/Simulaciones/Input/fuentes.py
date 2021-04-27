import numpy as np


def fuente_pde(x_grid, t_grid, Nx, Nt, sigma, dist, source):
    if source == 'gaussian':
        fuente_init = np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        fuente = [fuente_init] * Nt
        fuente = np.array(fuente)
        fuente = np.transpose(fuente)
        for j in range(Nx):
            fuente[j, :] = fuente[j, :]
        fuente = np.transpose(fuente)

    elif source == 'bigaussian':
        fuente_init_1 = np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))
        fuente_init_2 = -np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2))
        fuente_1 = [fuente_init_1] * Nt
        fuente_2 = [fuente_init_2] * Nt
        fuente_1 = np.array(fuente_1)
        fuente_1 = np.transpose(fuente_1)
        fuente_2 = np.array(fuente_2)
        fuente_2 = np.transpose(fuente_2)
        fuente = fuente_1 + fuente_2
        fuente = np.transpose(fuente)

    elif source == 'none':
        fuente = 0
    return fuente
