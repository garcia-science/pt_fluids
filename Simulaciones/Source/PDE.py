import numpy as np
from Simulaciones.Source.RK4 import *


def PDE_funcion(eq, PHI, i, j, dx, Nx_pasos, parametros, BC):
    if eq == 'Fisher':
        if BC == 'periodicas':
            if j == 0:
                F = (parametros[0] / dx ** 2) * (PHI[i, int(Nx_pasos) - 1] - 2 * PHI[i, j] + PHI[i, j + 1]) + parametros[1] * PHI[i, j] * (1 - PHI[i, j])
            elif j == int(Nx_pasos) - 1:
                F = (parametros[0] / dx ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, 0]) + parametros[1] * PHI[
                    i, j] * (1 - PHI[i, j])
            else:
                F = (parametros[0] / dx ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, j + 1]) + parametros[1] * PHI[
                    i, j] * (1 - PHI[i, j])
    if eq == 'KdV':
        if BC == 'periodicas':
            if j == 1:
                F = -(parametros[0] * PHI[i, j] * ((PHI[i, j + 1] - PHI[i, int(Nx_pasos) - 1]) / (2 * dx)) + parametros[1] * (
                            (-PHI[i, j - 2] + 2 * PHI[i, int(Nx_pasos) - 1] - 2 * PHI[i, j + 1] + PHI[i, j + 2]) / (8 * dx ** 3)))
            elif j == 0:
                F = -(parametros[0] * PHI[i, j] * ((PHI[i, j + 1] - PHI[i, int(Nx_pasos) - 1]) / (2 * dx)) + parametros[1] * (
                            (-PHI[i, int(Nx_pasos) - 2] + 2 * PHI[i, int(Nx_pasos) - 1] - 2 * PHI[i, j + 1] + PHI[i, j + 2]) / (8 * dx ** 3)))
            elif j == int(Nx_pasos) - 2:
                F = -(parametros[0] * PHI[i, j] * ((PHI[i, j + 1] - PHI[i, j - 1]) / (2 * dx)) + parametros[1] * (
                            (-PHI[i, j - 2] + 2 * PHI[i, j - 1] - 2 * PHI[i, j + 1] + PHI[i, 0]) / (8 * dx ** 3)))
            elif j == int(Nx_pasos) - 1:
                F = -(parametros[0] * PHI[i, j] * ((PHI[i, 0] - PHI[i, j - 1]) / (2 * dx)) + parametros[1] * (
                            (-PHI[i, j - 2] + 2 * PHI[i, j - 1] - 2 * PHI[i, 0] + PHI[i, 1]) / (8 * dx ** 3)))
            else:
                F = -(parametros[0] * PHI[i, j]*((PHI[i, j + 1] - PHI[i, j - 1])/(2*dx)) + parametros[1]*((-PHI[i, j - 2] + 2*PHI[i, j - 1] - 2*PHI[i, j + 1] + PHI[i, j + 2])/(8*dx**3)))
    if eq == 'Wave':
        if BC == 'periodicas':
            if j == 0:
                F = ((parametros[0] / dx) ** 2) * (PHI[i, int(Nx_pasos) - 1] - 2 * PHI[i, j] + PHI[i, j + 1])
            elif j == int(Nx_pasos) - 1:
                F = ((parametros[0] / dx) ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, 0])
            else:
                F = ((parametros[0] / dx) ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, j + 1])
    return F
