from procesos import *
from visualizacion import *
from matplotlib import pyplot as plt

if __name__ == '__main__':
    [X, T, PHI] = cargar_txt(experimental_data_file, data_file, X='X', T='T', PHI='PHI')
    visualizacion(X, T, PHI,
                  tipo="colormap", guardar='no',
                  path=experimental_data_file, file=data_file + '\datos_procesados',
                  nombre='plot')
    plt.show()