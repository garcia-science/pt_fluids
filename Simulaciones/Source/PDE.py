import numpy as np
from Simulaciones.Source.RK4 import *


def PDE_funciones(i , j, eq, campos, bordes, parametros, dx, Nx, Nt, x_grid, t_grid):
    if eq == 'pndls':
        U = campos[0]
        V = campos[1]
        if bordes == 'fixed':
            U_left = U_right = 0
            V_left = V_right = 0
            forzamiento = parametros[0]
            if j == 0:
                F = (U_left ** 2 + V_left ** 2) * U_left + parametros[2] * U_left + (parametros[3] * forzamiento[i, j] - parametros[4]) * V_left
                G = -((U_left ** 2 + V_left ** 2) * V_left + parametros[2] * V_left + (parametros[3] * forzamiento[i, j] + parametros[3]) * U_left)
            elif j == 1:
                F = parametros[1] * ((1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U_left) + (U[i, j] ** 2 + V[i, j] ** 2) \
                    * U[i, j] + parametros[2] * U[i, j] + (parametros[3] * forzamiento[i, j] - parametros[4]) * V[i, j])
                G = -parametros[1] * ((1 / dx ** 2) * (V[i, j + 1] - 2 * V[i, j] + V_left) + (U[i, j] ** 2 + V[i, j] ** 2)
                      * V[i, j] + parametros[2] * V[i, j] + (parametros[3] * forzamiento[i, j] + parametros[4]) * U[i, j])
            elif j == Nx - 1:
                F = (U_right ** 2 + V_right ** 2) * U_right + parametros[2] * U_right + (parametros[3] * forzamiento[i, j] - parametros[4]) * V_right
                G = -((U_right ** 2 + V_right ** 2) * V_right + parametros[2] * V_right + (parametros[3] * forzamiento[i, j] + parametros[4]) * U_right)
            elif j == Nx - 2:
                F = parametros[1] * ((1 / dx ** 2) * (U_right - 2 * U[i, j] + U[i, j - 1]) + (U[i, j] ** 2 + V[i, j] ** 2) \
                    * U[i, j] + parametros[2] * U[i, j] + (parametros[3] * forzamiento[i, j] - parametros[4]) * V[i, j])
                G = -parametros[1] * ((1 / dx ** 2) * (V_right - 2 * V[i, j] + V[i, j - 1]) + (U[i, j] ** 2 + V[i, j] ** 2)
                      * V[i, j] + parametros[2] * V[i, j] + (parametros[3] * forzamiento[i, j] + parametros[4]) * U[i, j])
            else:
                F = parametros[1] * ((1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1]) + (U[i, j] ** 2 + V[i, j] ** 2) \
                    * U[i, j] + parametros[2] * U[i, j] + (parametros[3] * forzamiento[i, j] - parametros[4]) * V[i, j])
                G = -parametros[1] * ((1 / dx ** 2) * (V[i, j + 1] - 2 * V[i, j] + V[i, j - 1]) + (U[i, j] ** 2 + V[i, j] ** 2)
                      * V[i, j] + parametros[2] * V[i, j] + (parametros[3] * forzamiento[i, j] + parametros[4]) * U[i, j])
            FUN = [F, G]
    elif eq == 'wave':
        U = campos[0]
        V = campos[1]
        if bordes == 'fixed':
            U_left = U_right = 0
            V_left = V_right = 0
            if j == 0:
                F = V_left
                G = 0
            elif j == 1:
                F = V[i, j]
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U_left)
            elif j == Nx - 1:
                F = V_right
                G = 0
            elif j == Nx - 2:
                F = V[i, j]
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U_right - 2 * U[i, j] + U[i, j - 1])
            else:
                F = -()
                G = (parametros[0] ** 2) * (1 / dx ** 2) * (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1])
            FUN = [F, G]
    elif eq == 'Fisher':
        if bordes== 'periodicas':
            U = campos[0]
            if j == 0:
                F = (parametros[0] / dx ** 2) * (U[i, Nx - 1] - 2 * U[i, j] + U[i, j + 1]) + parametros[1] * U[i, j] * (1 - U[i, j])
            elif j == Nx - 1:
                F = (parametros[0] / dx ** 2) * (U[i, j - 1] - 2 * U[i, j] + U[i, 0]) + parametros[1] * U[
                    i, j] * (1 - U[i, j])
            else:
                F = (parametros[0] / dx ** 2) * (U[i, j - 1] - 2 * U[i, j] + U[i, j + 1]) + parametros[1] * U[
                    i, j] * (1 - U[i, j])
            FUN = [F]
    elif eq == 'KdV':
        if bordes == 'periodicas':
            U = campos[0]
            if j == 1:
                F = -(parametros[0] * U[i, j] * ((U[i, j + 1] - U[i, Nx - 1]) / (2 * dx)) + parametros[1] * (
                            (-U[i, j - 2] + 2 * U[i, Nx - 1] - 2 * U[i, j + 1] + U[i, j + 2]) / (8 * dx ** 3)))
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
