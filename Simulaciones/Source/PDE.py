import numpy as np
from Simulaciones.Source.RK4 import *


def PDE_funciones(i , j, eq, campos, bordes, parametros, dx, Nx, Nt, x_grid, t_grid):
    if eq == 'pndls':
        U = campos[0]
        V = campos[1]
        if bordes == 'periodic':
            U_min = U[i, Nx - 1]
            U_max = U[i, 0]
            V_min = V[i, Nx - 1]
            V_max = V[i, 0]
            forzamiento = parametros[0]
            control = parametros[1]
            alpha = parametros[2]
            beta = parametros[3]
            gamma = parametros[4]
            mu = parametros[5]
            nu = parametros[6]
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
    elif eq == 'wave':
        U = campos[0]
        V = campos[1]
        if bordes == 'fixed':
            U_left = U_right = 0
            V_left = V_right = 0
            if j == 0:
                F = V_left
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U[i, j + 1] - 2 * U_left + (U_left - V_left * dx))
            elif j == 1:
                F = V[i, j]
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U_left)
            elif j == Nx - 1:
                F = V_right
                G = (parametros[0] ** 2) * (1 / dx ** 2) * ((U_right + V_right * dx) - 2 * U_right + U[i, j - 1])
            elif j == Nx - 2:
                F = V[i, j]
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U_right - 2 * U[i, j] + U[i, j - 1])
            else:
                F = V[i, j]
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1])
            FUN = [F, G]
        if bordes == 'periodic':
            U_min = U[i, Nx - 1]
            U_max = U[i, 0]
            if j == 0:
                F = V[i, j]
                G = ((parametros[0] / dx) ** 2) * (U_min - 2 * U[i, j] + U[i, j + 1])
            elif j == Nx - 1:
                F = V[i, j]
                G = ((parametros[0] / dx) ** 2) * (U[i, j - 1] - 2 * U[i, j] + U_max)
            else:
                F = V[i, j]
                G = ((parametros[0] / dx) ** 2) * (U[i, j - 1] - 2 * U[i, j] + U[i, j + 1])
            FUN = [F, G]
    elif eq == 'fisher':
        if bordes == 'periodic':
            U = campos[0]
            U_min = U[i, Nx - 1]
            U_max = U[i, 0]
            if j == 0:
                F = (parametros[0] / dx ** 2) * (U_min - 2 * U[i, j] + U[i, j + 1]) + parametros[1] * U[i, j] * (1 - U[i, j])
            elif j == Nx - 1:
                F = (parametros[0] / dx ** 2) * (U[i, j - 1] - 2 * U[i, j] + U_max) + parametros[1] * U[
                    i, j] * (1 - U[i, j])
            else:
                F = (parametros[0] / dx ** 2) * (U[i, j - 1] - 2 * U[i, j] + U[i, j + 1]) + parametros[1] * U[
                    i, j] * (1 - U[i, j])
            FUN = [F]
    elif eq == 'KdV':
        if bordes == 'periodic':
            U = campos[0]
            if j == 1:
                F = -(parametros[0] * U[i, j] * ((U[i, j + 1] - U[i, j - 1]) / (2 * dx)) + parametros[1] * (
                            (-U[i, Nx - 1] + 2 * U[i, j - 1] - 2 * U[i, j + 1] + U[i, j + 2]) / (8 * dx ** 3)))
            elif j == 0:
                F = -(parametros[0] * U[i, j] * ((U[i, j + 1] - U[i, Nx - 1]) / (2 * dx)) + parametros[1] * (
                            (-U[i, Nx - 2] + 2 * U[i, Nx - 1] - 2 * U[i, j + 1] + U[i, j + 2]) / (8 * dx ** 3)))
            elif j == Nx - 2:
                F = -(parametros[0] * U[i, j] * ((U[i, j + 1] - U[i, j - 1]) / (2 * dx)) + parametros[1] * (
                            (-U[i, j - 2] + 2 * U[i, j - 1] - 2 * U[i, j + 1] + U[i, 0]) / (8 * dx ** 3)))
            elif j == Nx - 1:
                F = -(parametros[0] * U[i, j] * ((U[i, 0] - U[i, j - 1]) / (2 * dx)) + parametros[1] * (
                            (-U[i, j - 2] + 2 * U[i, j - 1] - 2 * U[i, 0] + U[i, 1]) / (8 * dx ** 3)))
            else:
                F = -(parametros[0] * U[i, j]*((U[i, j + 1] - U[i, j - 1])/(2*dx)) + parametros[1]*(
                        (-U[i, j - 2] + 2 * U[i, j - 1] - 2 * U[i, j + 1] + U[i, j + 2])/(8*dx**3)))
            FUN = [F]
    return FUN
