import numpy as np


def variables(initial_condition):
    dt = float(input('Ingrese dt:'))
    t_min = int(input('Ingrese t_min:'))
    t_max = int(input('Ingrese t_max:'))
    print("Si acá aparece un error, debes colocar otro dt (cambiar esto por resto=0)")
    N = int((t_max-t_min)/dt)
    t = np.linspace(t_min, t_max, N)
    y = np.zeros(len(t))
    y[0] = initial_condition
    return dt, t, N, y