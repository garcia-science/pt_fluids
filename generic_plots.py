from procesos import *
from visualizacion import *
from numpy import genfromtxt

if __name__ == '__main__':
    ###   OPEN FILES   ###
    print('Abriendo archivos...')
    carpeta = select_directory(simulation_data_path)
    X = genfromtxt(carpeta + '\\' + 'X.txt')
    mod = genfromtxt(carpeta + '\\' + 'mod.txt')
    real = genfromtxt(carpeta + '\\' + 'real.txt')
    T = np.arange(len(real[:, 0]))

    plot = plt.pcolormesh(X, T, real, cmap='PiYG', shading='auto')
    cbar = plt.colorbar(plot, shrink=1)

    plt.savefig(carpeta + '\\' + 'real')
    plt.show()