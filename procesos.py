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
from tkinter import *
import tkinter as tk
from tkinter import filedialog
import os
from scipy.stats import linregress
from visualizacion import *


# NOMBRAR, GUARDAR Y CARGAR DATOS


def select_file(datos_path):
    root = tk.Tk()
    root.withdraw()
    carpeta = filedialog.askopenfilename(parent=root,
                                      initialdir=datos_path,
                                      title='Selecciones el archivo')
    return carpeta


def select_directory(datos_path):
    root = tk.Tk()
    root.withdraw()
    carpeta = filedialog.askdirectory(parent=root,
                                      initialdir=datos_path,
                                      title='Selecciones la carpeta')
    return carpeta


def crear_directorios_trabajo():
    root = tk.Tk()
    root.withdraw()

    def crear_directorio(path):
        if os.path.exists(path) == True:
            print('Este archivo ya existe')

        else:
            os.makedirs(path)
            print(path + ' creado')

    detection_parent_file = filedialog.askdirectory(parent=root,
                                                    initialdir='C:/',
                                                    title='Detección multiple')
    crear_directorio(detection_parent_file + '/mnustes_science/images/canned')
    crear_directorio(detection_parent_file + '/mnustes_science/images/img_lab')
    crear_directorio(detection_parent_file + '/mnustes_science/images/img_phantom')
    crear_directorio(detection_parent_file + '/mnustes_science/experimental_data')
    crear_directorio(detection_parent_file + '/mnustes_science/simulation_data')
    main_directory = detection_parent_file + '/mnustes_science'
    return main_directory


def guardar_txt(path, file, **kwargs): # upgradear a diccionario para nombre de variables
    if os.path.exists(path + file) == False:
        os.makedirs(path + file)
    for key, value in kwargs.items():
        np.savetxt(path + file + '\\' + key + ".txt", value)


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


def nombre_pndls_bigaussian(gamma, mu, nu, sigma, dist, fase):
    gamma_st = str(round(float(gamma), 3))
    mu_st = str(round(float(mu), 3))
    nu_st = str(round(float(nu), 3))
    sigma_st = str(round(float(sigma), 1))
    dist_st = str(round(float(dist), 1))
    fase_st = str(round(fase / np.pi, 2)) + 'pi'
    nombre = '\\bigaussian\\mu=' + mu_st + '\\gamma=' + gamma_st + '_nu=' + nu_st +'\\fase=' + fase_st +'\\sigma=' + sigma_st + '\\distancia=' + dist_st
    return nombre


# DETECCION

def canny_prueba(sigma):
    root = tk.Tk()
    root.withdraw()
    reference_image = filedialog.askopenfilename(parent=root,
                                                    initialdir="D:\mnustes_science",
                                                    title='Detección multiple')
    print(str(reference_image))
    im = cv2.imread(str(reference_image))
    REC = cv2.selectROI(im)
    rec = list(REC)
    imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
    imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
    canned = auto_canny(imBlur, sigma)
    cv2.imshow('Imagen de referencia', canned)
    cv2.waitKey(delay=0)
    cv2.destroyWindow('Imagen de referencia')


def canny_to_data():
    canned_path = 'D:\mnustes_science\images\canned'
    datos_path = 'D:\mnustes_science\experimental_data'
    root = tk.Tk()
    root.withdraw()
    detection_file = filedialog.askdirectory(parent=root,
                                                    initialdir=canned_path,
                                                    title='Selecciones la carpeta canny')
    if not detection_file:
        sys.exit('No se seleccionó ninguna carpeta')
    os.chdir(detection_file)
    parent_file_name = os.path.basename(detection_file)
    print('Se va a procesar la carpeta ' + detection_file)
    IMGs = os.listdir(canned_path + '\\single_file\\' + parent_file_name)
    X, T, PHI = datos_3d(IMGs, canned_path + '\\single_file\\' + parent_file_name, nivel='si')
    guardar_txt(datos_path, '\\single_file\\' + parent_file_name + '\\', X=X, T=T, PHI=PHI)


def deteccion_contornos(tipo, sigma, img_format, **kwargs):
    if tipo == 'multiple':
        root = tk.Tk()
        root.withdraw()
        detection_parent_file = filedialog.askdirectory(parent=root,
                                                        initialdir="D:\mnustes_science",
                                                        title='Detección multiple')
        if not detection_parent_file:
            sys.exit('No se seleccionó ningún archivo')
        os.chdir(detection_parent_file)
        detection_files = os.listdir()
        parent_file_name = os.path.basename(detection_parent_file)
        print('Se va a procesar la carpeta ' + str(parent_file_name))
        canned_path = 'D:\mnustes_science\images\canned'
        datos_path = 'D:\mnustes_science\experimental_data'

        reference_image = filedialog.askopenfilename(parent=root,
                                                        initialdir=detection_files,
                                                        title='Seleccionar imagen de referencia')
        recs = ROI_select(reference_image)
        for name in detection_files:
            print('Procesando ' + str(name) + ' (' + str(detection_files.index(name)) + '/' + str(len(detection_files)) + ')')
            if img_format == 'jpg':
                deteccion_jpg(detection_parent_file + '\\' + name, canned_path + '\\' + parent_file_name + '\\' + name,
                              recs, sigma)
            elif img_format == 'tiff':
                deteccion_tiff(detection_parent_file + '\\' + name, canned_path + '\\' + parent_file_name + '\\' + name,
                              recs, sigma)
            IMGs = os.listdir(canned_path + '\\' + parent_file_name + '\\' + name)
            X, T, PHI = datos_3d(IMGs, canned_path + '\\' + parent_file_name + '\\' + name)
            guardar_txt(datos_path, '\\' + parent_file_name + '\\' + name, X=X, T=T, PHI=PHI)
    elif tipo == 'single_file':
        root = tk.Tk()
        root.withdraw()
        zero_file = filedialog.askdirectory(parent=root,
                                                        initialdir="D:\mnustes_science",
                                                        title='Seleccione la carpeta del cero')
        if not zero_file:
            sys.exit('No se seleccionó ningún archivo')

        detection_file = filedialog.askdirectory(parent=root,
                                                 initialdir="D:\mnustes_science",
                                                 title='Seleccione la carpeta para detección')
        if not detection_file:
            sys.exit('No se seleccionó ningún archivo')
        os.chdir(detection_file)
        parent_file_name = os.path.basename(detection_file)
        os.chdir(detection_file)
        zero_name = os.path.basename(zero_file)
        canned_path = 'D:\mnustes_science\images\canned'
        datos_path = 'D:\mnustes_science\experimental_data'
        print('Se va a procesar la carpeta ' + detection_file)
        reference_image = filedialog.askopenfilename(parent=root,
                                                     initialdir=detection_file,
                                                     title='Seleccionar imagen de referencia')
        recs = ROI_select(reference_image)
        if 'file_name' not in kwargs:
            file_name = 'default'
        else:
            file_name = kwargs['file_name']

        if img_format == 'jpg':
            deteccion_jpg(zero_file, canned_path + '\\' + file_name + '\\' + parent_file_name + '\\' + zero_name, recs, sigma)
        elif img_format == 'tiff':
            deteccion_tiff(zero_file, canned_path + '\\' + file_name + '\\' + parent_file_name + '\\' + zero_name, recs, sigma)
        IMGs = os.listdir(canned_path + '\\' + file_name + '\\' + parent_file_name + '\\' + zero_name)
        X, T, ZERO = datos_3d(IMGs, canned_path + '\\' + file_name + '\\' + parent_file_name + '\\' + zero_name)
        guardar_txt(datos_path, '\\' + file_name + '\\' + parent_file_name, ZERO=ZERO)
        if img_format == 'jpg':
            deteccion_jpg(detection_file, canned_path + '\\' + file_name + '\\' + parent_file_name, recs, sigma)
        elif img_format == 'tiff':
            deteccion_tiff(detection_file, canned_path + '\\' + file_name + '\\' + parent_file_name, recs, sigma)
        IMGs = os.listdir(canned_path + '\\' + file_name + '\\' + parent_file_name)
        X, T, PHI = datos_3d(IMGs, canned_path + '\\' + file_name + '\\' + parent_file_name)
        guardar_txt(datos_path, '\\' + file_name + '\\' + parent_file_name , X=X, T=T, PHI=PHI)
    return X, T, PHI


def auto_canny(image, sigma):
    if sigma == 'fixed':
        lower = 100
        upper = 200
    else:
        v = np.median(image)
        lower = int(max(0, (1.0 - sigma) * v))
        upper = int(min(255, (1.0 + sigma) * v))
    edged = cv2.Canny(image, lower, upper)
    return edged


def deteccion_jpg(file_i, file_o, REC, sigma):
    IMGs = os.listdir(file_i)  # lista de nombres de archivos en la carpeta indicada
    im = cv2.imread(file_i + '/cam000000.jpg')
    rec = list(REC)
    imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
    imBlur = cv2.GaussianBlur(imCrop, (7, 7), 0)
    ddepth = cv2.CV_16S
    scale = 1
    delta = 0
    grad_x = cv2.Sobel(imBlur, ddepth, 1, 0, ksize=3, scale=scale, delta=delta, borderType = cv2.BORDER_DEFAULT)
    grad_y = cv2.Sobel(imBlur, ddepth, 0, 1, ksize=3, scale=scale, delta=delta, borderType = cv2.BORDER_DEFAULT)
    abs_grad_x = cv2.convertScaleAbs(grad_x)
    abs_grad_y = cv2.convertScaleAbs(grad_y)
    grad = cv2.addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0)
    edges = auto_canny(grad, sigma)
    if os.path.exists(file_o) == True:
        print('Este archivo de CANNY ya existe, ¿desea eliminarlo y continuar? (y/n)')
        a = str(input())
        if a == 'y':
            shutil.rmtree(file_o)
        elif a == 'n':
            sys.exit("Proceso terminado, cambie de carpeta")

    os.makedirs(file_o)
    cv2.imwrite(os.path.join(file_o, IMGs[0]), edges)
    for i in range(1, len(IMGs)):
        im = cv2.imread(file_i + '\\' + IMGs[i])
        imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
        imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
        # edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur, sigma)
        cv2.imwrite(os.path.join(file_o, IMGs[i]), edges)
    return IMGs


def deteccion_tiff(file_i, file_o, REC, sigma):
    IMGs = os.listdir(file_i)  # lista de nombres de archivos en la carpeta indicada
    im = cv2.imread(file_i + '/cam000000.tif')
    rec = list(REC)
    imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
    imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
    edges = auto_canny(imBlur, sigma)
    if os.path.exists(file_o) == True:
        print('Este archivo de CANNY ya existe, ¿desea eliminarlo y continuar? (y/n)')
        a = str(input())
        if a == 'y':
            shutil.rmtree(file_o)
        elif a == 'n':
            sys.exit("Proceso terminado, cambie de carpeta")

    os.makedirs(file_o)
    cv2.imwrite(os.path.join(file_o, IMGs[0]), edges)
    for i in range(1, len(IMGs)):
        im = cv2.imread(file_i + '\\' + IMGs[i])
        imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
        imBlur = cv2.GaussianBlur(imCrop, (3, 3), 0)
        # edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur, sigma)
        cv2.imwrite(os.path.join(file_o, IMGs[i]), edges)
    return IMGs


def ROI_select(path):
    im = cv2.imread(path)
    RECs = cv2.selectROI(im)
    return RECs


# IMAGENES A DATOS
def phi_t(IMGs, file_o, l):
    img = cv2.imread(file_o + '\\' + IMGs[l], 0)
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
    return phi, cols


def datos_3d(IMGS, FILE_OUT):
    PHI = []
    T = []
    N_imgs = len(IMGS)
    if N_imgs == 1:
        phi, cols = phi_t(IMGS, FILE_OUT, 0)
        t = [0]
        PHI.append(phi)
        T.append(t)
    else:
        for i in range(1, N_imgs):
            phi, cols = phi_t(IMGS, FILE_OUT, i)
            t = [i]
            PHI.append(phi)
            T.append(t)
    X = np.arange(1, cols)
    Y = np.array(T)
    Z = np.array(PHI)
    return X, Y, Z


# PROCESOS DE DATOS

def drift_velocity(T_per, X_mm, Z_mm, window_l, window_u, t_inicial, t_final):
    ###   DEFINIENDO COSAS, VENTANA INICIAL E INTERVALO TEMPORAL A ANALIZAR  ###
    L_wind = window_u - window_l

    ###   ENCONTRANDO MAXIMOS   ###
    t_array = []
    x_array = []
    for i in range(t_inicial, t_final):
        j = window_l + np.argmax(Z_mm[i, window_l:window_u])
        t_array.append(T_per[i])
        x_array.append(X_mm[j])
        window_l = int(j - L_wind / 2)
        window_u = int(j + L_wind / 2)
    t_np = np.array(t_array)
    x_np = np.array(x_array)

    ###   REGRESIÓN LINEAL   ###
    linear_fit = linregress(t_array, x_array)
    x_fit = linear_fit.slope * t_np + linear_fit.intercept
    return t_np, x_np, x_fit, linear_fit


def zero_fix(z_limit, mode, cargar, *args):
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    if mode == 'zero':
        if cargar == 'si':
            [X, T, PHI, zero] = cargar_txt(carpeta, '', X='X', T='T', PHI='PHI', ZERO='ZERO')
        elif cargar == 'no':
            zero = cargar_txt(carpeta, '', ZERO='ZERO')
            [X, T, PHI] = [args[0], args[1], args[2]]
        ZERO = np.ones((len(PHI[:, 0]), len(PHI[0, :])))
        for i in range(len(T)):
            ZERO[i, :] = zero
        Z = PHI - ZERO
        Z = np.array(Z)
        guardar_txt(carpeta, '', Z=Z)
        visualizacion(X, T, Z, tipo='colormap', guardar='si', path=carpeta,
                      file='', nombre='espaciotiempo_mean', cmap='seismic', vmin=-z_limit, vzero=0, vmax=z_limit)
        plt.close()
    elif mode == 'mean':
        if cargar == 'si':
            [X, T, PHI] = cargar_txt(carpeta, '', X='X', T='T', PHI='PHI')
        elif cargar == 'no':
            [X, T, PHI] = [args[0], args[1], args[2]]
        Z = nivel_mean(PHI, X, T)
        Z = np.array(Z)
        guardar_txt(carpeta, '', Z=Z)
        visualizacion(X, T, Z, tipo='colormap', guardar='si', path=carpeta,
                      file='', nombre='espaciotiempo_filt', cmap='seismic', vmin=-z_limit, vzero=0, vmax=z_limit)
        plt.close()
    return carpeta, X, T, Z

def nivel_mean(PHI, X, T):
    mean = np.mean(PHI[:, 0])
    #PHI = filtro_superficie(PHI, 3, 'X')
    MEAN = mean * np.ones((len(PHI[:, 0]), len(PHI[0, :])))
    Z = PHI - MEAN
    mmin = Z[0, 0]
    mmax = Z[0, -1]
    pend = mmax - mmin
    nivels = []
    for i in range(len(X)):
        y_i = (pend / len(X)) * X[i]
        nivels.append(y_i)
    nivels = np.array(nivels)
    Z_new = []
    for i in range(len(T)):
        Z_new_i = Z[i, :] - nivels
        Z_new_i = Z_new_i.tolist()
        Z_new.append(Z_new_i)
    return Z_new


def field_envelopes(X, T, Z, carpeta):
    def envelopes(s):
        q_u = np.zeros(s.shape)
        q_l = np.zeros(s.shape)
        u_x = [0, ]
        u_y = [s[0], ]
        l_x = [0, ]
        l_y = [s[0], ]
        for k in range(1, len(s) - 1):
            if (np.sign(s[k] - s[k - 1]) == 1) and (np.sign(s[k] - s[k + 1]) == 1):
                u_x.append(k)
                u_y.append(s[k])
            if (np.sign(s[k] - s[k - 1]) == -1) and ((np.sign(s[k] - s[k + 1])) == -1):
                l_x.append(k)
                l_y.append(s[k])
        u_x.append(len(s) - 1)
        u_y.append(s[-1])
        l_x.append(len(s) - 1)
        l_y.append(s[-1])
        u_p = interp1d(u_x, u_y, kind='linear', bounds_error=False, fill_value=0.0)
        l_p = interp1d(l_x, l_y, kind='linear', bounds_error=False, fill_value=0.0)
        for k in range(0, len(s)):
            q_u[k] = u_p(k)
            q_l[k] = l_p(k)
        q_u = q_u.tolist()
        q_l = q_l.tolist()
        return q_u, q_l
    A = np.zeros((len(T), len(X)))
    B = np.zeros((len(T), len(X)))
    for i in range(len(X)):
        print(i)
        s = Z[:, i]
        q_u, q_l =envelopes(s)
        A[:, i] = q_u
        B[:, i] = q_l
    guardar_txt(carpeta, '', A=A, B=B)
    visualizacion(X, T, A, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='A_plot', cmap='seismic')
    plt.close()
    visualizacion(X, T, B, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='B_plot', cmap='seismic')
    plt.close()


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


def sin(x, A, w, phase):
    fun = A * np.sin(w * x + phase)
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


def fit_sin(X, Y):
    X = np.array(X)
    Y = np.array(Y)
    popt, pcov = curve_fit(sin, X, Y)
    fit = sin(X, *popt)
    return fit, popt

def ajuste_altura(X, T, Z, threshold, n):
    Z_ajustado = np.ones((len(Z[:, 0]), len(Z[0, :])))
    for i in range(len(T)):
        for j in range(len(X)):
            if Z[i, j] > threshold:
                Z_ajustado[i, j] = threshold + (Z[i, j] - threshold) * n
            else:
                Z_ajustado[i, j] = Z[i, j]
    return Z_ajustado