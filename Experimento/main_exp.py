from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *
import cv2
import os
from scipy.optimize import curve_fit

if __name__ == '__main__':
    ############ DETECTOR DE CONTORNOS ###########

    #recs = ROI_select(reference_file)
    #deteccion(path, file_in, file_out, recs, sigma = 0.37)

    ########## GUARDADO Y CARGA DE DATOS #########
    #IMGs = os.listdir(path + file_out)
    #X, T, PHI = datos_3d(IMGs, path, file_out, 100, filtro = 'no')
    #datos = [X, T, PHI]
    #guardar_txt(datos, data_path, data_file)
    [X, T, PHI] = cargar_txt(3, data_path, data_file)
    X = np.array([i - X[-1] / 2 for i in X])

    ########## PROCESOS Y OTROS #########
    PHI_filtrado = filtro_superficie(PHI, 10, 'YX')
    mean, std = proyeccion_desvesta(PHI_filtrado)
    PHI_filtrado = nivel(PHI_filtrado, mean)
    PHI_proy, frec, power_density = proyeccion_maximos(PHI_filtrado)
    std = (PHI_proy[np.argmax(PHI_proy)] / std[np.argmax(std)]) * std
    envelope, puntos_x, puntos_y = envelope(X, std, 'linear')
    fit, popt = fit_gauss_sin(X, std)
    #fit, popt = gaussian_fit(X, std, X)
    #################### CONVERSION A CM ##################
    arrays = [X, PHI_proy, std, fit]
    [X, PHI_proy, std, fit] = resize_arrays(5, arrays)

    #################### PLOTEO ##################
    #fig1 = plot_ZXT(X, Y, Z, guardar='no', nombre='datos_3D_01', titulo='Datos 3D')
    #fig2 = color_map(X, Y, Z, guardar='no', nombre, titulo, xlabel=r'$x$ (Espacio)', ylabel=r'$\psi(x)$ (Altura)')
    #plot_XY(X, envelope, guardar='no', nombre='', titulo='Gráfico interesante', xlabel=r'$x$ (Espacio)', ylabel=r'$\psi(x)$ (Altura)')

    multiple_XY(X, [std, fit], guardar='no', nombre='comparacion_01', titulo='Estandar v/s Fourier',
               xlabel=r'$x$ (Espacio)', ylabel=r'$\psi(x)$ (Altura)')
    #plt.scatter(X, std)
    #fig4 = DOS(frec, power_density, "log", "Periodograma",guardar='no', nombre='periodograma_01', titulo='Periodograma', xlabel=r'$\omega$ (Frecuencias)', ylabel=r'Intensidad')
    plt.show()



