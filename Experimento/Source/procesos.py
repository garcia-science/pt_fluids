import cv2
from scipy.signal import lfilter, filtfilt
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

def filtro_array(n, funcion):
    n = 10  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    phi_filtered = filtfilt(b, a, funcion)
    return phi_filtered

def phi_t(main_path, IMGs , file_o, filtro, l, m, nivel_0, nivel_f):
    img = cv2.imread(main_path + file_o + '\\' + IMGs[l], 0)
    rows, cols = img.shape
    # 1) Va por columna viendo que pixel es blanco,
    # lo guarda y dividiendo por el número de
    # pixeles en una columna para luego anotar
    # su valor entre 0 y 1 en un array.
    # 2) El contador cuenta las filas al revés (upsi),
    # al final se invierte el array para que quede derecho.
    phi = []
    i = cols - 1
    while i != 0:
        j = rows - 1
        while j != 0:
            n = 0
            k = img[j, i]
            if k == 255:
                phi_i = 1 - (j / rows)
                phi.append(phi_i)
                j = 0
                n = 1
            elif k != 255:
                j = j - 1
            if j == 1 and n == 0:
                if not phi:
                    phi_i = 0.45
                    print("tengo un hoyo al principio")
                    phi.append(phi_i)
                    j = j - 1
                else:
                    phi_i = phi[-1]
                    phi.append(phi_i)
                    j = j - 1
        i = i - 1
    x = []
    for i in range(cols - 1):
        x.append(i)
    phi.reverse()
    if filtro == 'si':
        phi_F = filtro_array(m, phi)
    elif filtro == 'no':
        phi_F = np.array(phi)
    if nivel_0 == 0:
        nivel_0 = phi_F[0]
        nivel_f = phi_F[-1]
    phi_F = phi_F - [nivel_0] * (cols - 1)  # Centrando en el cero
    pend = (nivel_f - nivel_0)/cols #Corrector de nivel
    recta = []
    for a in x:
        recta.append(a * pend)
    phi_F = phi_F - np.array(recta)
    return phi_F, cols, nivel_0, nivel_f

def datos_3d(IMGS, PATH, FILE_OUT, m, filtro):
    PHI = []
    T = []
    nivel_0 = 0
    nivel_f = 0
    for i in range(1, len(IMGS)):
        print(i)
        phi, cols, nivel_0, nivel_f = phi_t(PATH, IMGS, FILE_OUT, filtro, i, m, nivel_0, nivel_f)
        t = [i]
        PHI.append(phi)
        T.append(t)
    X = np.arange(1, cols)
    Y = np.array(T)
    Z = np.array(PHI)
    return PHI, T, X, Y, Z

def proyeccion(PHI):
    PHIT = PHI.transpose()
    rows, cols = PHIT.shape
    PHIT_proy = np.zeros(rows)
    for i in range(cols-1):
        PHIT_proy = PHIT_proy + np.absolute(PHIT[:, i])
    PHIT_proy = (1/cols)*PHIT_proy
    return PHIT_proy

def proyeccion_maximos(Z):
    phi_max = np.argmax(Z[0, :])
    frecuencias, power_density = signal.periodogram(Z[:, phi_max])
    max_element = np.argmax(power_density)
    periodo = 1 / frecuencias[max_element]
    print("El periodo es "+str(periodo))
    max_int = np.argmax(Z[0, 0:int(periodo)])
    max_int = int(max_int)
    A = np.array([Z[max_int, :], Z[max_int, :]])
    for i in range(len(Z[:, phi_max])):
        if int(max_int * (2 + i/2)) < len(Z[:, 0]):
            A_i = np.array(np.absolute([Z[int(max_int * (2 + i/2)), :]]))
            A = np.append(A, A_i, axis=0)
    A = A[1:-1]
    PHIT_proy = proyeccion(A)
    return PHIT_proy, frecuencias, power_density

def proyeccion_desvesta(Z):
    N_x = len(Z[0, :])
    N_t = len(Z[:, 0])
    print(N_x, N_t)
    mean = np.array([])
    std = np.array([])
    for i in range(N_x):
            mean_i = np.array([np.mean(Z[:, i])])
            mean = np.append(mean, mean_i)
            std_i = np.array([np.std(Z[:, i])])
            std = np.append(std, std_i) #falta factr multiplicativo para esto
    return mean, std

def envelope(X, Y, tipo):
    u_x = [X[0]]
    u_y = [Y[0]]
    for i in range(1, len(Y)-1):
        if np.sign(Y[i]-Y[i-1]) == 1 and np.sign(Y[i]-Y[i+1]) == 1:
            u_x.append(X[i])
            u_y.append(Y[i])
    u_x.append(X[-1])
    u_y.append(Y[-1])
    print(u_x)
    print(u_y)
    u_p = interp1d(u_x, u_y, kind = tipo,bounds_error = False, fill_value=0.0)
    #https://towardsdatascience.com/basic-curve-fitting-of-scientific-data-with-python-9592244a2509
    q_u = []
    for k in range(0, len(Y)):
        q_u.append(u_p(k))
    return q_u, u_p

def resize_array(cm_px, array):
    array_cm = [i * cm_px for i in array]
    return array_cm






















