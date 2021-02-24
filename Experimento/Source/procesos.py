import cv2
from scipy.signal import lfilter, filtfilt
import numpy as np


def filtro_array(n, funcion):
    n = 100  # the larger n is, the smoother curve will be
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
    phi_F = phi_F - [(nivel_0)] * (cols - 1)  # Centrando en el cero
    pend = (nivel_f - nivel_0)/cols #Corrector de nivel
    recta = []
    for a in x:
        recta.append(a * pend)
    phi_F = phi_F - np.array(recta)
    return phi_F, cols, nivel_0, nivel_f

def datos_3d(IMGS, PATH, FILE_OUT, filtro):
    PHI = []
    T = []
    nivel_0 = 0
    nivel_f = 0
    for i in range(1, len(IMGS)):
        print(i)
        phi, cols, nivel_0, nivel_f = phi_t(PATH, IMGS, FILE_OUT, filtro, i, 10, nivel_0, nivel_f)
        t = [i]
        PHI.append(phi)
        T.append(t)
    X = np.arange(1, cols)
    Y = np.array(T)
    Z = np.array(PHI)
    return PHI, T, X, Y, Z

def proyeccion(PHI): #cuenta de más
    PHIT = PHI.transpose()
    rows, cols = PHIT.shape
    PHIT_proy = np.zeros(rows)
    for i in range(cols-1):
        PHIT_proy = PHIT_proy + np.absolute(PHIT[:, i])
    PHIT_proy = (1/cols)*PHIT_proy
    return PHIT_proy

#def periodograma()