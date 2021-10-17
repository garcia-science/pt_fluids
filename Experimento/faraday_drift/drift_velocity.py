from procesos import *
from visualizacion import *
from numpy import genfromtxt


if __name__ == '__main__':
    ###   OPEN FILES   ###
    datos_path = 'D:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    T_per = genfromtxt(carpeta + '\\' + 'T_per.csv', delimiter=',')
    X_mm = genfromtxt(carpeta + '\\' + 'X_mm.csv', delimiter=',')
    Z_mm = genfromtxt(carpeta + '\\' + 'Z_mm.csv', delimiter=',')

    ###   DEFINIENDO COSAS, VENTANA INICIAL E INTERVALO TEMPORAL A ANALIZAR  ###
    T = np.arange(len(T_per))
    X = np.arange(len(X_mm))

    window_l = 524
    window_u = 576
    t_inicial = 10000
    t_final = 30000
    L_wind = window_u - window_l

    ###   ENCONTRANDO MAXIMOS   ###
    i_array = []
    j_array = []
    tx_array = []
    x_array = []
    for i in range(t_inicial, t_final):
        j = window_l + np.argmax(Z_mm[i, window_l:window_u])
        tx_array.append(T_per[i])
        x_array.append(X_mm[j])
        window_l = int(j - L_wind / 2)
        window_u = int(j + L_wind / 2)
    tx_np = np.array(tx_array)
    x_np = np.array(x_array)

    ###   REGRESIÓN LINEAL   ###
    linear_fit = linregress(tx_array, x_array)
    x_fit = linear_fit.slope * tx_np + linear_fit.intercept

    ###   PLOTS   ###
    plt.plot(x_fit, tx_array, linewidth=3, c='g')
    plot = plt.pcolormesh(X_mm, T_per, Z_mm, cmap='Reds', shading='auto')
    cbar = plt.colorbar(plot, shrink=1)
    plt.show()