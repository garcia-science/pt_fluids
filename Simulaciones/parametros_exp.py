import numpy as np
if __name__ == '__main__':
    unidad = 0.1
    g = 9790 * unidad
    l_x = 160 * unidad
    l_y = 16 * unidad
    d = 20 * unidad
    n = 5
    m = 1
    a_ang = 7
    f = 15.00

    a_mm = (28 + a_ang / 12) * unidad
    l_iny = l_x / 3
    sigma = l_iny / 2.3482
    k_x = np.pi * n / l_x
    k_y = np.pi * m / l_y
    k = np.sqrt(k_x ** 2 + k_y ** 2)
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau) / (2 * np.pi)
    w = f / 2
    GAMMA = 4 * w ** 2 * a_mm
    print("w = " + str(w))
    print("w_1 = " + str(w_1))
    print("k_x = " + str(k_x))
    print("k_y = " + str(k_y))
    print("k = " + str(k))
    print("tau = " + str(tau))

    alpha = (1 / (4 * k ** 2)) * (1 + k * d * ((1 - tau ** 2) / tau)) #término difusivo
    beta = (k ** 2 / 64) * (6 * tau ** 2 - 5 + 16 * tau ** (-2) - 9 * tau ** (-4)) #término no lineal
    gamma = GAMMA / (4 * g) #amplitud adimensional
    nu = 0.5 * ((w / w_1) ** 2 - 1)
    mu = 0.0125
    print("ADIMENSIONALES")
    print("alpha = " + str(alpha) + "\nbeta = " + str(beta) + "\ngamma = " + str(gamma) + "\nnu = " + str(nu))