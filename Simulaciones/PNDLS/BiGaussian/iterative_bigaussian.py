from procesos import *
from visualizacion import *
from Simulaciones.Input.inicializacion import *
from Simulaciones.Recursos.evolucion import *

if __name__ == '__main__':
    alpha = 1
    beta = 1
    gamma = 0.28
    mu = 0.1
    nu = 0.32
    sigma_forcing_1 = sigma_forcing_2 = sigma_forcing = 3
    distancia = 20
    fase = np.pi
    L = 200
    dx = 0.5
    dt = 0.001
    T_final = 600
    iteration_1 = ['distancia', distancia, 21, 1]
    iteration_2 = ['sigma_forcing', sigma_forcing, 4, 0.1]
    iterative_bigaussian_adimensional(iteration_1, iteration_2, alpha, beta, gamma, mu, nu, sigma_forcing_1, sigma_forcing_2, distancia, fase, L, dx, dt, T_final)