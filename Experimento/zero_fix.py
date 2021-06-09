import numpy as np
from scipy.interpolate import interp1d
from procesos import *
from visualizacion import *


if __name__ == '__main__':
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    Z = zero_fix(carpeta, 20, 'filt')
