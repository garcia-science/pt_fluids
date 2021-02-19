from Source.ODE import *

def RK4(F, y, dt, t, eq, parametros):
    for i in range(0, len(t)-1):
        k1 = F
        k2 = F + 0.5*dt*k1
        k3 = F + 0.5*dt*k2
        k4 = F + dt*k3
        y[i+1] = y[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6
        F = funcion(eq, y, parametros)
    return y