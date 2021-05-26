from procesos import *
from visualizacion import *
from scipy.interpolate import interp1d

if __name__ == '__main__':
    datos_path = 'F:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    [zero, X, T, PHI] = cargar_txt(carpeta, '', ZERO='ZERO', X='X', T='T', PHI='PHI')
    Z = nivel_mean(PHI, X, T)
    #ZERO = np.ones((len(PHI[:, 0]), len(PHI[0, :])))
    #for i in range(len(T)):
    #    ZERO[i, :] = zero
    #Z = PHI - ZERO
    Z = np.array(Z)
    guardar_txt(carpeta, '', Z=Z)
    visualizacion(X, T[0:200], Z[0:200, :], tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='diagrama_espaciotiempo', cmap='viridis')