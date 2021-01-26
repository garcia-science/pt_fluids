import cv2
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import lfilter, filtfilt
from PIL import Image

path = r'C:\Users\mnustes_science\PT_fluids'
file_in = '\IMG'
file_out = '\IMG_canned'
reference_file = r'C:\Users\mnustes_science\PT_fluids\IMG'

def auto_canny(image, sigma):
    v = np.median(image)
    lower = int(max(0, (1.0 - sigma) * v))
    upper = int(min(255, (1.0 + sigma) * v))
    edged = cv2.Canny(image, lower, upper)
    return edged


def deteccion(main_path, file_i, file_o, REC, sigma):
    IMGs = os.listdir(main_path + file_i)  # lista de nombres de archivos en la carpeta indicada
    if __name__ == '__main__':
        im = cv2.imread(main_path + file_i + '\cam0001.jpg')
        rec = list(REC)
        imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
        imBlur = cv2.GaussianBlur(imCrop, (9, 9), 6, 6, cv2.BORDER_DEFAULT)
        # edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur, sigma)
        cv2.imwrite(os.path.join(main_path + file_o, IMGs[0]), edges)
        for i in range(1, len(IMGs)):
            im = cv2.imread(main_path + file_i + '\\' + IMGs[i])
            imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
            imBlur = cv2.GaussianBlur(imCrop, (5, 5), 6, 6, cv2.BORDER_DEFAULT)
            # edges = cv2.Canny(imBlur,10,200)
            edges = auto_canny(imBlur, sigma)
            cv2.imwrite(os.path.join(main_path + file_o, IMGs[i]), edges)
    return IMGs


def ROI_select(path):
    if __name__ == '__main__':
        im = cv2.imread(path + '\cam0001.jpg')
        fromCenter = False
        REC = cv2.selectROI(im)
    return REC


def filtro_array(n, funcion):
    n = 100  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    phi_filtered = filtfilt(b, a, funcion)
    return phi_filtered


def phi_t(main_path, file_o, filtro, l, m):
    if __name__ == '__main__':
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
        phi_filtered = filtro_array(m, phi)
    elif filtro == 'no':
        phi_filtered = phi
    return phi_filtered, cols

#rec = ROI_select(reference_file)
#deteccion(path, file_in, file_out, rec, 4)

IMGs = os.listdir(path + file_out)
PHI = []
T = []
for i in range(1, len(IMGs)):
    print(i)
    phi, cols = phi_t(path, file_out, 'si', i, 10)
    t = [i]
    PHI.append(phi)
    T.append(t)
np.array(PHI)
X = np.arange(1, cols)
Y = np.array(T)
X, Y = np.meshgrid(X, Y)
Z = np.array(PHI)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
ax.set_zlim(0, 1)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()