import numpy as np


def PDE_funciones(i, j, eq,  campos, bordes, dx, Nx, **kwargs):
    if eq == 'pndls' or 'pndls_exp':
        U = campos[0]
        V = campos[1]
        if bordes == 'periodic':
            U_min = U[i, Nx - 1]
            U_max = U[i, 0]
            V_min = V[i, Nx - 1]
            V_max = V[i, 0]
            forzamiento = kwargs['forzamiento']
            control = kwargs['control']
            mu = kwargs['mu']
            alpha = kwargs['alpha']
            beta = kwargs['beta']
            nu = kwargs['nu']
            gamma = kwargs['gamma']
            if j == 0:
                F = control * ((gamma * forzamiento[i, j] - mu) * U[i, j] + nu * V[i, j] + (alpha * (1 / dx) ** 2) * (
                            V_min - 2 * V[i, j] + V[i, j + 1]) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * V[i, j])
                G = - control * ((gamma * forzamiento[i, j] + mu) * V[i, j] + nu * U[i, j] + (alpha * (1 / dx) ** 2) * (
                            U_min - 2 * U[i, j] + U[i, j + 1]) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * U[i, j])
            elif j == Nx - 1:
                F = control * ((gamma * forzamiento[i, j] - mu) * U[i, j] + nu * V[i, j] + (alpha * (1 / dx) ** 2) * (
                            V[i, j - 1] - 2 * V[i, j] + V_max) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * V[i, j])
                G = - control * ((gamma * forzamiento[i, j] + mu) * V[i, j] + nu * U[i, j] + (alpha * (1 / dx) ** 2) * (
                            U[i, j - 1] - 2 * U[i, j] + U_max) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * U[i, j])
            else:
                F = control * ((gamma * forzamiento[i, j] - mu) * U[i, j] + nu * V[i, j] + (alpha * (1 / dx) ** 2) * (
                            V[i, j - 1] - 2 * V[i, j] + V[i, j + 1]) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * V[i, j])
                G = - control * ((gamma * forzamiento[i, j] + mu) * V[i, j] + nu * U[i, j] + (alpha * (1 / dx) ** 2) * (
                            U[i, j - 1] - 2 * U[i, j] + U[i, j + 1]) + beta * (U[i, j] ** 2 + V[i, j] ** 2) * U[i, j])
            FUN = [F, G]
    #print('alpha=' + str(alpha) + ' beta=' + str(beta) + ' mu=' + str(mu) + ' gamma=' + str(gamma) + ' nu=' + str(
    #           nu) + ' control=' + str(control))
    return FUN

def funcion(eq, y, parametros):
    if eq == "logistica":
        F = parametros[0] * y * (1-y/parametros[1])
    return F
