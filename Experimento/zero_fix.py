import numpy as np
from scipy.interpolate import interp1d
from procesos import *
from visualizacion import *


if __name__ == '__main__':
    carpeta, X, T, Z = zero_fix(20, 'mean', 'si')
    visualizacion(X, T, Z, tipo='colormap', guardar='si', path=carpeta,
                 file='', nombre='faraday', cmap='seismic', vmin=-20, vzero=0, vmax=20)
    plt.show()