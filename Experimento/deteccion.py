from directorios import *
from procesos import *
from visualizacion import *
import winsound

if __name__ == '__main__':
    sigma = 1.2
    #=crear_directorios_trabajo()
    #canny_prueba(sigma)
    deteccion_contornos('single_file', sigma, 'jpg')
    winsound.PlaySound('C:\\Users\\rariv\\Downloads\\alarm.wav', winsound.SND_FILENAME)