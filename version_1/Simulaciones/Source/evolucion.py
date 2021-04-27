import numpy as np
from Simulaciones.Source.RK4 import *
from Simulaciones.Source.PDE import *

def PDE_temporal(PHI_inicial, dx, dt, Nt_pasos, Nx_pasos, eq, parametros, BC):
    if eq == 'Fisher' or eq == 'KdV':
        PHI = RK4_PDE(PHI_inicial, dx, dt, Nt_pasos, Nx_pasos, eq, parametros, BC)
    elif eq =="Wave":
        for i in range(0, int(Nt_pasos) - 1):
            print('dt = ' + str(i))
            for j in range(0, int(Nx_pasos) - 1):  # poner if de condiciones de borde
                F = PDE_funcion(eq, PHI_inicial, i, j, dx, Nx_pasos, parametros, BC)
                k1 = F
                k2 = F + 0.5 * dt * k1
                k3 = F + 0.5 * dt * k2
                k4 = F + dt * k3
                PHI_inicial[i + 1, j] = PHI_inicial[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return PHI