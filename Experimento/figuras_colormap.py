from procesos import *
from visualizacion import *

if __name__ == '__main__':
    [X, T, PHI] = cargar_txt(experimental_data_file, '\\07-05-2021\\f=15.22_a=6.0', X='X', T='T', PHI='PHI')
    mean = np.mean(PHI[:, 0])
    #PHI = filtro_superficie(PHI, 20, 'X')
    MEAN = mean * np.ones((len(PHI[:, 0]), len(PHI[0, :])))
    Z = PHI - MEAN
    mmin = Z[0, 0]
    mmax = Z[0, -1]
    pend = mmax - mmin
    nivels = []
    for i in range(len(X)):
        y_i = (pend / len(X)) * X[i]
        nivels.append(y_i)
    nivels = np.array(nivels)
    Z_new = []
    for i in range(len(T)):
        Z_new_i = Z[i, :] - nivels
        Z_new_i = Z_new_i.tolist()
        Z_new.append(Z_new_i)
    Z_new = np.array(Z_new)
    print(len(T))
    #D = []
    #n=100
    #for i in range(len(T)):
    #    if i*n < 12000:
    #        aaa = np.argmax(Z[i * n, :])
    #        bbb = np.argmax(Z[i * n:(i + 1) * n, aaa])
    #        print(aaa)
    #        print(bbb)
    #        C = Z[bbb, :]*n
    #        D.append(C)
    #    else:
    #        break
    #D = np.array(D)
    #visualizacion(X, np.arange(0, len(D[:, 0]), 1), D, tipo='colormap', guardar='si', path=experimental_data_file, file='\\07-05-2021\\f=15.22_a=9.4', nombre='plot', cmap='seismic')
    #

    #visualizacion(X, [Z_new[0, :]], tipo="2D_multiple", guardar='no', path=experimental_data_file,
    #              file='\\07-05-2021\\f=15.22_a=6.0', nombre='plot')
    #visualizacion(T[0:100], a[0:100, :], tipo='2D', guardar='no', path=experimental_data_file,
    #              file='\\07-05-2021\\f=15.22_a=9.4', nombre='plot')
    visualizacion(X, T, Z_new, tipo='colormap', guardar='si', path=experimental_data_file,
                  file='\\07-05-2021\\f=15.22_a=6.0', nombre='plot', cmap='seismic')
    plt.show()