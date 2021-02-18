import numpy as np

def funcion(eq, y, parametros):
    if eq == "logistica":
        F = parametros[0]*y*(1-y/parametros[1])
    elif eq == "ola":
        F = np.exp(-y)
    return F