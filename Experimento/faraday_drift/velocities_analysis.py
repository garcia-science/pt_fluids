from procesos import *
from visualizacion import *
from numpy import genfromtxt


if __name__ == '__main__':
    ampd = 11.8

    ###   OPEN FILES   ###
    print('Abriendo archivos...')
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    dvelocity_input = genfromtxt(carpeta + '\\' + '/inputs/dvelocity_input_' + str(ampd) + '.csv', delimiter=',')
    info_dvelocity = genfromtxt(carpeta + '\\' + 'info_dvelocity_' + str(ampd) + '.csv', delimiter=',')
    mean = np.mean(info_dvelocity[:, 2])
    std = np.std(info_dvelocity[:, 2])

    x = [0, ampd]
    y = [0, mean]
    y_err = [0, std]

    plt.errorbar(x, y, yerr=y_err, marker='o', ls='', capsize=5, capthick=1, ecolor='black', color='red')
    plt.show()