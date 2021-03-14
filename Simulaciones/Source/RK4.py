from Simulaciones.Source.ODE import *
from Simulaciones.Source.PDE import *

def RK4_ODE(F, y, dt, t, eq, parametros):
    for i in range(0, len(t)-1):
        print(i)
        k1 = F
        k2 = F + 0.5 * dt * k1
        k3 = F + 0.5 * dt * k2
        k4 = F + dt * k3
        y[i+1] = y[i] + dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])/6
        F = funcion(eq, y, parametros)
    return y


def RK4_PDE(campos, bordes, parametros, dx, dt, Nx, Nt, eq, x_grid, t_grid):
    for i in range(Nt - 1):
        print('dt = ' + str(i))
        for j in range(Nx ):
            FUN = PDE_funciones(i, j, eq, campos, bordes, parametros, dx, Nx, Nt, x_grid, t_grid)
            K = [0, 0, 0, 0]
            for k in range(len(FUN)):
                K[0] = FUN[k]
                K[1] = FUN[k] + 0.5 * dt * K[0]
                K[2] = FUN[k] + 0.5 * dt * K[1]
                K[3] = FUN[k] + dt * K[2]
                campos[k][i + 1, j] = campos[k][i, j] + dt * (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6
    return campos