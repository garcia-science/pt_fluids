from procesos import *
from visualizacion import *
from numpy import genfromtxt
import os

if __name__ == '__main__':
    ampd = 11.8

    ###   OPEN FILES   ###
    print('Abriendo archivos...')
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    T_per = genfromtxt(carpeta + '\\' + 'T_per.csv', delimiter=',')
    X_mm = genfromtxt(carpeta + '\\' + 'X_mm.csv', delimiter=',')
    Z_mm = genfromtxt(carpeta + '\\' + 'Z_mm.csv', delimiter=',')
    dvelocity_input = [[230, 260, 10000, 25000], [332, 372, 10000, 25000], [420, 460, 10000, 25000], [540, 570, 10000, 25000], [620, 660, 10000, 25000]]

    print('Comenzando iteración')
    info_dvelocity = []
    fits = []
    for i in range(len(dvelocity_input)):
        t_np, x_np, x_fit, linear_fit = drift_velocity(T_per, X_mm, Z_mm, dvelocity_input[i][0], dvelocity_input[i][1], dvelocity_input[i][2], dvelocity_input[i][3])
        fits_i = [t_np, x_fit]
        info_dvelocity_i = [ampd, i, linear_fit.slope, linear_fit.intercept]
        info_dvelocity.append(info_dvelocity_i)
        fits.append(fits_i)
        print('Drift velocity ' + str(i + 1) + '/' + str(len(dvelocity_input)))
    dvelocity_input_np = np.array(dvelocity_input)
    info_dvelocity_np = np.array(info_dvelocity)
    os.chdir(carpeta)
    os.chdir('../')
    cwd = os.getcwd()
    dvelocities_file = cwd + '/dvelocities_info'
    inputs_file = dvelocities_file + '/inputs'
    if not os.path.exists(dvelocities_file):
        os.makedirs(dvelocities_file)
    if not os.path.exists(inputs_file):
        os.makedirs(inputs_file)
    np.savetxt(inputs_file + '\\' + 'dvelocity_input_' + str(ampd) + '.csv', info_dvelocity, delimiter=',')
    np.savetxt(dvelocities_file + '\\' + 'info_dvelocity_' + str(ampd) + '.csv', info_dvelocity, delimiter=',')

    plt.plot(fits[0][1], fits[0][0], linewidth=3, c='g')
    plt.plot(fits[1][1], fits[1][0], linewidth=3, c='g')
    plt.plot(fits[2][1], fits[2][0], linewidth=3, c='g')
    plt.plot(fits[3][1], fits[3][0], linewidth=3, c='g')
    plt.plot(fits[4][1], fits[4][0], linewidth=3, c='g')

    plot = plt.pcolormesh(X_mm, T_per, Z_mm, cmap='Reds', shading='auto')
    cbar = plt.colorbar(plot, shrink=1)
    plt.savefig(carpeta + '\\' + 'drift_velocities')
    plt.show()