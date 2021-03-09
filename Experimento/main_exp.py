from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *
import cv2
import os

if __name__ == '__main__':
    ############ DETECTOR DE CONTORNOS ###########
    #recs = ROI_select(reference_file)
    cm_px = escala(reference_file, 4)


    #deteccion(path, file_in, file_out, recs, sigma = 0.37)

    ############## IMAGENES A NUMPY ##############
    #IMGs = os.listdir(path + file_out)
    #PHI, T, X, Y, Z = datos_3d(IMGs, path, file_out, 100, filtro = 'si')

    ########## GUARDADO Y CARGA DE DATOS #########
    #datos = [X, Y, Z, PHI]
    #guardar_txt(datos)

    [X, Y, Z, PHI] = cargar_txt(4)
    FILT = np.zeros((len(Z[:, 0]), len(Z[0, :])))
    for i in range(len(Z[0, :])):
        filtered = filtro_array(10, Z[:, i])
        FILT[:, i] = filtered




    ########### PROYECCION DE MAXIMOS ############
    PHIT_proy, frec, power_density = proyeccion_maximos(FILT)


    ##########ENSAYOS###########

    mean, std = proyeccion_desvesta(FILT)
    envelope, puntos = envelope(X, std, 'linear')


    X_cm = np.array(resize_array(cm_px, X))
    envelope_cm = np.array(resize_array(cm_px, envelope))
    std_cm = np.array(resize_array(cm_px, std))
    PHIT_proy_cm = np.array(resize_array(cm_px, PHIT_proy))
    mean_cm = np.array(resize_array(cm_px, mean))
    E = envelope_cm - 1 * mean_cm
    S = std_cm - 1 * mean_cm
    PH = PHIT_proy_cm - 1 * mean_cm
    '''ARREGLAR CORRECCION DE NIVEL, PROBLEMA CON ALTURA, Y OTROS'''
    #################### PLOTEO ##################
    #fig1 = plot_ZXT(X, Y, Z)
    #fig2 = color_map(X, Y, FILT)
    plot_XY(X_cm, E)
    plot_XY(X_cm, S)
    plot_XY(X_cm, PH)
    #plot_XY(X_cm, mean_cm)
    #fig4 = DOS(frec, power_density, "log", "Periodograma")
    plt.show()
