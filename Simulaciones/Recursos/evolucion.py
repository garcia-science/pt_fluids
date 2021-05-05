from Simulaciones.Recursos.ecuaciones import *


def RK4_ODE(F, y, dt, t, eq, parametros):
    for i in range(0, len(t)-1):
        print(i)
        k1 = F
        k2 = F + 0.5 * dt * k1
        k3 = F + 0.5 * dt * k2
        k4 = F + dt * k3
        y[i+1] = y[i] + dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])/6
        F = funcion(eq, y, parametros)
    return y


def RK4_PDE(eq, campos, bordes, dx, dt, Nx, Nt, **kwargs):
    forzamiento = kwargs['forzamiento']
    control = kwargs['control']
    mu = kwargs['mu']
    if 'n' and 'forcing_amp' and 'forcing_freq' and 'sigma' and 'profundidad' not in kwargs:
        alpha = 1
        beta = 1
        nu = kwargs['nu']
        gamma = kwargs['gamma']
    elif 'nu' and 'gamma' and 'alpha' and 'beta' not in kwargs:
        g = 9806.65
        n = kwargs['n']
        a = kwargs['forcing_amp']
        f = kwargs['forcing_freq']
        w = 2 * f
        l_inyeccion = 2 * kwargs['sigma']
        d = kwargs['profundidad']
        k = 2 * np.pi * n / l_inyeccion
        tau = np.tanh(k * d)
        GAMMA = 4 * w ** 2 * a
        w_1 = np.sqrt(g * k * tau)
        alpha = (1 / (4 * k ** 2)) * (1 + k * d * ((1 - tau ** 2) / tau))
        beta = (k ** 2 / 64) * (6 * tau ** 2 - 5 + 16 ** tau ** (-2) - 9 * tau ** (-4))
        gamma = GAMMA / (4 * g)
        nu = 0.5 * ((w / w_1) ** 2 - 1)
    for i in range(Nt - 1):
        if i % int(1 / dt) == 0:
            print(str(round((i / (Nt - 1)), 3) * 100)+' %')
        for j in range(Nx):
            FUN = PDE_funciones(i, j, eq, campos, bordes, dx, Nx, control=control, alpha=alpha, beta=beta, gamma=gamma, mu=mu, nu=nu, forzamiento=forzamiento)
            K = [0, 0, 0, 0]
            for k in range(len(FUN)):
                K[0] = FUN[k]
                K[1] = FUN[k] + 0.5 * dt * K[0]
                K[2] = FUN[k] + 0.5 * dt * K[1]
                K[3] = FUN[k] + dt * K[2]
                campos[k][i + 1, j] = campos[k][i, j] + dt * (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6
    return campos