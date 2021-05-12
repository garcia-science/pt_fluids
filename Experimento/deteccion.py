from directorios import *
from procesos import *
from visualizacion import *


if __name__ == '__main__':
    sigma = 3
    #canny_prueba(sigma)
    deteccion_contornos('multiple', sigma)