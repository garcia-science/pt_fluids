import numpy as np


def fuente_pde(x_grid, t_grid, Nx, Nt, source):
    if source == 'gaussian':
        sigma = 20
        frec = 0.75
        fuente_init = np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        fuente = [fuente_init] * Nt
        fuente = np.array(fuente)
        fuente = np.transpose(fuente)
        for j in range(Nx):
            fuente[j, :] = fuente[j, :] * np.sin(frec * t_grid)
        fuente = np.transpose(fuente)

    elif source == 'bigaussian':
        sigma = 6
        frec = 0.5
        dist = 40
        fuente_init_1 = np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))
        fuente_init_2 = np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2))
        fuente_1 = [fuente_init_1] * Nt
        fuente_2 = [fuente_init_2] * Nt
        fuente_1 = np.array(fuente_1)
        fuente_1 = np.transpose(fuente_1)
        fuente_2 = np.array(fuente_2)
        fuente_2 = np.transpose(fuente_2)
        for j in range(Nx):
            fuente_1[j, :] = fuente_1[j, :] * np.sin(frec * t_grid)
            fuente_2[j, :] = fuente_2[j, :] * np.sin(frec * (t_grid + np.pi))
        fuente = fuente_1 + fuente_2
        fuente = np.transpose(fuente)

    elif source == 'none':
        fuente = 0
    return fuente
