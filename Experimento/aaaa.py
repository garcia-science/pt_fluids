import numpy as np
from scipy.interpolate import interp1d
from procesos import *
from visualizacion import *

if __name__ == '__main__':
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    #Z = zero_fix(carpeta, 20, 'mean')
    [X, T, Z] = cargar_txt(carpeta, '', X='X', T='T', Z='Z')
    Z_ajustado = ajuste_altura(X, T, Z, 6, 0.4)
    #visualizacion(X, T, Z_ajustado, tipo='colormap', guardar='no', path=carpeta,
    #              file='', nombre='strobo_plot_2_filt2', cmap='seismic')
    guardar_txt(carpeta, '', X=X, Z_ajustado=Z_ajustado)
    #plt.plot(Z[115, :])
    #plt.plot(Z_ajustado[115, :])
    #plt.show()
