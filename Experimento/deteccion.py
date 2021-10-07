from directorios import *
from procesos import *
from visualizacion import *
import winsound

if __name__ == '__main__':
    sigma = 'fixed'
    project_file = 'faraday_drift_02'
    X, T, PHI = deteccion_contornos('single_file', sigma, 'jpg', file_name=project_file)