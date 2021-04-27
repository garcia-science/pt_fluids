import cv2
from scipy.signal import filtfilt
import numpy as np
import os
import shutil
from scipy import signal
import sys
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from directorios import *
from visualizacion import *
from scipy.optimize import curve_fit


# NOMBRAR, GUARDAR Y CARGAR DATOS


def guardar_txt(path, file, **kwargs): # upgradear a diccionario para nombre de variables
    for key, value in kwargs.items():
        np.savetxt(path + file + '\\' + key +".txt", value)


def cargar_txt(path, file, **kwargs): # upgradear a diccionario para nombre de variables
    array = []
    for key, values in kwargs.items():
        array_i = np.loadtxt(path + file + '\\' + key + ".txt")
        array.append(array_i)
    return array


def nombre_pndls_estandar(**kwargs):
    mu = kwargs['mu']
    if 'n' and 'forcing_amp' and 'forcing_freq' and 'profundidad' not in kwargs:
        nu = kwargs['nu']
        gamma = kwargs['gamma']
        mu_st = str(round(float(mu), 3))
        gamma_st = str(round(float(gamma), 3))
        nu_st = str(round(float(nu), 3))
        mu_splited = mu_st.split('.')
        gamma_splited = gamma_st.split('.')
        nu_splited = nu_st.split('.')
        mu_name = mu_splited[0] + mu_splited[1]
        gamma_name = gamma_splited[0] + gamma_splited[1]
        nu_name = nu_splited[0] + nu_splited[1]
        nombre = '\\gaussian\mu=' + mu_name + '\gamma=' + gamma_name + '_nu=' + nu_name
    elif 'alpha' and 'beta' and 'nu' and 'gamma' not in kwargs:
        d = kwargs['profundidad']
        n = kwargs['n']
        a = kwargs['forcing_amp']
        w = kwargs['forcing_freq']
        d_name = str(round(float(d), 2))
        n_name = str(n)
        a_name = str(round(float(a), 2))
        w_name = str(round(float(w), 2))
        nombre = '\\gaussian_exp\\d=' + d_name + '\\n=' + n_name + '\\f=' + w_name + '_a=' + a_name
    return nombre


def nombre_pndls_bigaussian(gamma, nu, sigma, dist, fase):
    gamma_st = str(round(float(gamma), 3))
    nu_st = str(round((float(nu), 3)))
    sigma_st = str(round(float(sigma), 1))
    dist_st = str(round(float(dist), 1))
    gamma_splited = gamma_st.split('.')
    nu_splited = nu_st.split('.')
    sigma_splited = sigma_st.split('.')
    dist_splited = dist_st.split('.')
    gamma_name = gamma_splited[0] + gamma_splited[1]
    nu_name = nu_splited[0] + nu_splited[1]
    sigma_name = sigma_splited[0] + sigma_splited[1]
    dist_name = dist_splited[0] + dist_splited[1]
    nombre = '\\bigaussian\gamma=' + gamma_name + '_nu=' + nu_name +'\\fase=' + fase +'\\sigma=' + sigma_name + '_distancia=' + dist_name
    return nombre


# DETECCION


def auto_canny(image, sigma):
    v = np.median(image)
    lower = int(max(0, (1.0 - sigma) * v))
    upper = int(min(255, (1.0 + sigma) * v))
    edged = cv2.Canny(image, lower, upper)
    return edged


def deteccion(main_path, file_i, file_o, REC, sigma):
    IMGs = os.listdir(main_path + file_i)  # lista de nombres de archivos en la carpeta indicada
    im = cv2.imread(main_path + file_i + '\cam0001.jpg')
    rec = list(REC)
    imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
    imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
    edges = auto_canny(imBlur, sigma)
    if os.path.exists(main_path + file_o) == True:
        print('Este archivo de CANNY ya existe, ¿desea eliminarlo y continuar? (y/n)')
        a = str(input())
        if a == 'y':
            shutil.rmtree(main_path + file_o)
        elif a == 'n':
            sys.exit("Proceso terminado, cambie de carpeta")

    os.makedirs(main_path + file_o)
    cv2.imwrite(os.path.join(main_path + file_o, IMGs[0]), edges)
    for i in range(1, len(IMGs)):
        im = cv2.imread(main_path + file_i + '\\' + IMGs[i])
        imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
        imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
        # edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur, sigma)
        cv2.imwrite(os.path.join(main_path + file_o, IMGs[i]), edges)
    return IMGs


def ROI_select(path):
    im = cv2.imread(path + '\\cam0001.jpg')
    fromCenter = False
    RECs = cv2.selectROI(im)
    return RECs


# IMAGENES A DATOS
def phi_t(main_path, IMGs, file_o, l, nivel):
    img = cv2.imread(main_path + file_o + '\\' + IMGs[l], 0)
    rows, cols = img.shape
    phi = []
    i = cols - 1
    while i != 0:
        j = rows - 1
        while j != 0:
            n = 0
            k = img[j, i]
            if k == 255:
                phi_i = rows - j
                phi.append(phi_i)
                j = 0
                n = 1
            elif k != 255:
                j = j - 1
            if j == 1 and n == 0:
                if not phi:
                    phi_i = 0.5 * rows
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
    if nivel == 'si':
        nivel = phi[0]
    phi = [i - nivel for i in phi]  # Centrando en el cero
    return phi, cols, nivel


def datos_3d(IMGS, PATH, FILE_OUT, nivel):
    PHI = []
    T = []
    N_imgs = len(IMGS)
    for i in range(1, N_imgs):
        print(str(round((i / N_imgs) * 100, 1)) + '% procesado')
        phi, cols, nivel = phi_t(PATH, IMGS, FILE_OUT, i, nivel)
        t = [i]
        PHI.append(phi)
        T.append(t)
    X = np.arange(1, cols)
    Y = np.array(T)
    Z = np.array(PHI)
    return X, Y, Z


# PROCESOS DE DATOS


def filtro_array(n, funcion):
    # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    phi_filtered = filtfilt(b, a, funcion)
    return phi_filtered


def filtro_superficie(Z, intensidad, sentido):
    X_len = len(Z[:, 0])
    Y_len = len(Z[0, :])
    FILT = np.zeros((X_len, Y_len))
    if sentido == 'X':
        for i in range(X_len):
            filtered = filtro_array(intensidad, Z[i, :])
            FILT[i, :] = filtered
    elif sentido == 'Y':
        for i in range(Y_len):
            filtered = filtro_array(intensidad, Z[:, i])
            FILT[:, i] = filtered
    elif sentido == 'XY':
        for i in range(X_len):
            filtered = filtro_array(intensidad, Z[i, :])
            FILT[i, :] = filtered
        for i in range(Y_len):
            filtered = filtro_array(intensidad, FILT[:, i])
            FILT[:, i] = filtered
    elif sentido == 'YX':
        for i in range(Y_len):
            filtered = filtro_array(intensidad, Z[:, i])
            FILT[:, i] = filtered
        for i in range(X_len):
            filtered = filtro_array(intensidad, FILT[i, :])
            FILT[i, :] = filtered
    return FILT


def proyeccion_maximos(Z):
    def proyeccion(PHI):
        PHIT = PHI.transpose()
        rows, cols = PHIT.shape
        PHIT_proy = np.zeros(rows)
        for i in range(cols - 1):
            PHIT_proy = PHIT_proy + np.absolute(PHIT[:, i])
        PHIT_proy = (1 / cols) * PHIT_proy
        return PHIT_proy
    phi_inicial = Z[0, :]
    phi_max = np.argmax(phi_inicial) # tomar el argumento del máximo valor del primer contorno
    maximo_temporal = Z[:, phi_max] # array de como se comporta el maximo encontrado en el tiempo
    frecuencias, power_density = signal.periodogram(maximo_temporal) # periodograma del array anterior
    max_element = np.argmax(power_density) # toma la frecuencia que corresponde al ajuste sinosidal
    periodo = 1 / frecuencias[max_element] # periodo asociado a la frecuencia
    max_int = np.argmax(Z[0:int(periodo), phi_max]) # encuentra el máximo en el primer periodo del
    max_int = int(max_int)
    A = []
    for i in range(1, 2 * int(len(Z[:, phi_max])/periodo)):
        if int(max_int * i/2) < len(Z[:, 0]):
            A_i = np.absolute(Z[int(max_int * i / 2), :])
            A.append(A_i)
    A = A[1:]
    A_np = np.array(A)
    PHIT_proy = proyeccion(A_np)
    return PHIT_proy, frecuencias, power_density


def proyeccion_desvesta(Z):
    N_x = len(Z[0, :])
    mean = np.array([])
    std = np.array([])
    for i in range(N_x):
            mean_i = np.array([np.mean(Z[:, i])])
            mean = np.append(mean, mean_i)
            std_i = np.array([np.std(Z[:, i])])
            std = np.append(std, std_i)
    std = std - std[0]
    return mean, std


def resize_arrays_referenced(cm, arrays, file):
    arrays_cm = []
    cm_px = escala(file, cm)
    for j in range(len(arrays)):
        array_cm = [i * cm_px for i in arrays[j]]
        arrays_cm.append(array_cm)
    return arrays_cm


def resize_arrays(cm, dx, arrays):
    arrays_cm = []
    cm_px = cm / dx
    for j in range(len(arrays)):
        array_cm = [i * cm_px for i in arrays[j]]
        arrays_cm.append(array_cm)
    return arrays_cm


def campos_ligeros(campos, n, Nt, Nx, T):
    t_ligero = np.linspace(0, T, int(Nt / n))
    campos_light = []
    for k in range(len(campos)):
        campo_ligero = np.zeros((int(Nt / n), Nx))
        for i in range(0, len(campos[k][:, 0]) - 1, n):
            campo_ligero[int(i / n), :] = campos[k][i, :]
        campos_light.append(campo_ligero)
    return campos_light, t_ligero


def escala(path, cm):
    dim = ROI_select(path)
    DIM = list(dim)
    cm_px = cm/DIM[2]
    print(cm_px)
    return cm_px


def nivel(Z, mean):
    for i in range(len(Z[:, 0])):
        Z[i, :] = Z[i, :] - mean
    return Z

# AJUSTES


def gauss_sin(x, a, sigma, L):
    fun = a * np.exp(- 0.5 * ((x / sigma) ** 2)) * (np.sin(2 * np.pi * x / L)) ** 2
    return fun


def gauss(x, a, sigma):
    fun = a * np.exp(- 0.5 * ((x / sigma) ** 2))
    return fun


def fit_gauss_sin(X, Y):
    X = np.array(X)
    Y = np.array(Y)
    a_min = np.amax(Y)
    sigma_max = 4 * np.amax(X)
    L_max = 0.25 * sigma_max
    popt, pcov = curve_fit(gauss_sin, X, Y, bounds=([a_min, 0, 0], [np.inf, sigma_max, L_max]))
    fit = gauss_sin(X, *popt)
    return fit, popt

















