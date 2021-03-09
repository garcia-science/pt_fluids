from Simulaciones.Source.ODE import *
from Simulaciones.Source.PDE import *

def RK4_ODE(F, y, dt, t, eq, parametros):
    for i in range(0, len(t)-1):
        print(i)
        k1 = F
        k2 = F + 0.5*dt*k1
        k3 = F + 0.5*dt*k2
        k4 = F + dt*k3
        y[i+1] = y[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6
        F = funcion(eq, y, parametros)
    return y
def RK4_PDE(PHI, dx, dt, Nt_pasos, Nx_pasos, eq, parametros, BC):
    for i in range(0, int(Nt_pasos)-1):
        #print('dt = '+str(i))
        for j in range(0, int(Nx_pasos)-1):#poner if de condiciones de borde
            F = PDE_funcion(eq, PHI, i, j, dx, Nx_pasos, parametros, BC)
            k1 = F
            k2 = F + 0.5 * dt * k1
            k3 = F + 0.5 * dt * k2
            k4 = F + dt * k3
            PHI[i + 1, j] = PHI[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return PHI