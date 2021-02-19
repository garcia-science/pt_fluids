import numpy as np
import matplotlib.pyplot as plt
from Experimento.Visualizacion.tres_dimensiones import plot_ZXT
if __name__ == '__main__':
    #Parametros de grilla
    dx = 0.1
    dt = 0.1
    x_min = -10
    x_max = 10
    T = 200
    Nx_pasos = (x_max - x_min)/dx
    Nt_pasos = T/dt
    print('tenemos '+str(Nx_pasos)+' pasos espaciales y '+str(Nt_pasos)+' pasos temporales')

    #Parametros ecuacion de fisher
    D = 0.01
    r = 0.1

    #Preparacion de grilla y condicion inicial
    x_grid = np.linspace(x_min, x_max, int(Nx_pasos))
    t_grid = np.linspace(0, T, int(Nt_pasos))
    phi_inicial = np.exp(-10*x_grid**2)
    #plt.plot(x_grid,phi_inicial)
    #plt.show()

    #Ecuacion
    PHI = np.zeros((int(Nx_pasos), int(Nt_pasos)-1))
    PHI = np.matrix.transpose(PHI)
    PHI = np.append([phi_inicial], PHI, axis=0)
    print(PHI.shape)
    print(phi_inicial.shape)
    print(t_grid.shape)
    #XX, TT = np.meshgrid(x_grid, t_grid)
    #fig = plot_ZXT(XX, TT, PHI)
    #plt.show()
    #print(PHI)
    for i in range(0, int(Nt_pasos)-1):
        print(i)
        for j in range(0, int(Nx_pasos)-1):#poner if de condiciones de borde
            if j == 0:
                F = (D/dx**2)*(PHI[i, int(Nx_pasos)-1] - 2*PHI[i, j] + PHI[i, j+1]) + r*PHI[i, j]*(1-PHI[i, j])
                k1 = F
                k2 = F + 0.5 * dt * k1
                k3 = F + 0.5 * dt * k2
                k4 = F + dt * k3
                PHI[i+1, j] = PHI[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            elif j == int(Nx_pasos)-1:
                F = (D / dx ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, 0]) + r * PHI[i, j] * (1 - PHI[i, j])
                k1 = F
                k2 = F + 0.5 * dt * k1
                k3 = F + 0.5 * dt * k2
                k4 = F + dt * k3
                PHI[i + 1, j] = PHI[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            else:
                F = (D / dx ** 2) * (PHI[i, j - 1] - 2 * PHI[i, j] + PHI[i, j + 1]) + r * PHI[i, j] * (1 - PHI[i, j])
                k1 = F
                k2 = F + 0.5 * dt * k1
                k3 = F + 0.5 * dt * k2
                k4 = F + dt * k3
                PHI[i + 1, j] = PHI[i, j] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    #print(PHI)
    XX, TT = np.meshgrid(x_grid, t_grid)
    fig = plot_ZXT(XX, TT, PHI)
    #fig = plt.imshow(PHI, interpolation='bilinear')
    plt.show()
    #for j in range(0,Nt_pasos):
    #    for k in range(0, Nx_pasos):





