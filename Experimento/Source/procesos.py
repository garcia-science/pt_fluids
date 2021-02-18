import cv2
from scipy.signal import lfilter, filtfilt
import numpy as np


def filtro_array(n, funcion):
    n = 100  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    phi_filtered = filtfilt(b, a, funcion)
    return phi_filtered


def phi_t(main_path, IMGs , file_o, filtro, l, m):
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
        phi_F = phi
    return phi_F, cols

def datos_3d(IMGS, PATH, FILE_OUT, filtro):
    PHI = []
    T = []
    for i in range(1, len(IMGS)):
        print(i)
        phi, cols = phi_t(PATH, IMGS, FILE_OUT, filtro, i, 10)
        t = [i]
        PHI.append(phi)
        T.append(t)
        X = np.arange(1, cols)
        Y = np.array(T)
        X_grid, Y_grid = np.meshgrid(X, Y)
        Z_np = np.array(PHI)
    return PHI, T, X_grid, Y_grid, Z_np