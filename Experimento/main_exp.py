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
    #guardar_txt(datos)
    [X, T, PHI] = cargar_txt(3)

    ########## PROCESOS Y OTROS #########
    PHI_filtrado = filtro_superficie(PHI, 5, 'YX')
    mean, std = proyeccion_desvesta(PHI_filtrado)
    PHI_filtrado = nivel(PHI_filtrado, mean)
    PHI_proy, frec, power_density = proyeccion_maximos(PHI_filtrado)
    std = (PHI_proy[np.argmax(PHI_proy)] / std[np.argmax(std)]) * std
    envelope, puntos_x, puntos_y = envelope(X, std, 'linear')



    #################### CONVERSION A CM ##################
    arrays = [X, PHI_proy, std, envelope]
    [X, PHI_proy, std, envelope] = resize_arrays(4, arrays)


    #################### CONVERSION A CM ##################
    def func(x, a, b, c):
        return a * np.exp(-(x - b) ** 2 / (2 * c) ** 2)
    plt.plot(X, PHI_proy, label='data')
    popt, pcov = curve_fit(func, X, envelope)
    plt.plot(X, func(X, *popt), 'r-',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()

    #################### PLOTEO ##################
    #fig1 = plot_ZXT(X, Y, Z)
    #fig2 = color_map(X, T, PHI_filtrado)
    plot_XY(X, PHI_proy)
    plot_XY(X, std)
    plot_XY(X, envelope)
    #fig4 = DOS(frec, power_density, "log", "Periodograma")
    plt.show()



