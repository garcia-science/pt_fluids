from Input.parametros import parametros
from Input.variables import *
from Source.ODE import *
from Source.RK4 import *
from matplotlib import pyplot as plt

if __name__ == '__main__':
    eqn, params, y_0 = parametros()
    h, T, n, x = variables(y_0)
    f = funcion(eqn, x, params)
    y_t = RK4_ODE(f, x, h, T, eqn, params)
    plt.plot(T, y_t)
    plt.show()
