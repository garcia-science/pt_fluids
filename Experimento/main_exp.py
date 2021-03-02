from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *

if __name__ == '__main__':
    ############ DETECTOR DE CONTORNOS ###########
    #recs = ROI_select(reference_file)
    #deteccion(path, file_in, file_out, recs, sigma = 0.37)

    ############## IMAGENES A NUMPY ##############
    #IMGs = os.listdir(path + file_out)
    #PHI, T, X, Y, Z = datos_3d(IMGs, path, file_out, filtro = 'si')

    ########## GUARDADO Y CARGA DE DATOS #########
    #datos = [X, Y, Z, PHI]
    #guardar_txt(datos)
    [X, Y, Z, PHI] = cargar_txt(4)

    ########### PROYECCION DE MAXIMOS ############
    PHIT_proy, frec, power_density = proyeccion_maximos(Z)

    #################### PLOTEO ##################
    #fig1 = plot_ZXT(X, Y, Z)
    #fig2 = color_map(X, Y, Z)
    fig3 = plot_XY(X, PHIT_proy)
    #fig4 = DOS(frec, power_density, "log", "Periodograma")
    plt.show()