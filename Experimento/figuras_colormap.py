from procesos import *
from visualizacion import *
from scipy.interpolate import interp1d

if __name__ == '__main__':
    datos_path = 'D:\mnustes_science\experimental_data'
    root = tk.Tk()
    root.withdraw()
    carpeta = filedialog.askdirectory(parent=root,
                                             initialdir=datos_path,
                                             title='Selecciones la carpeta a procesar')
    [X, T, PHI] = cargar_txt(carpeta, '', X='X', T='T', PHI='PHI')
    mean = np.mean(PHI[:, 0])
    #PHI = filtro_superficie(PHI, 3, 'X')
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

    Z = np.array(Z_new)
    guardar_txt(carpeta, '', Z=Z)
    field_envelopes(X, T, Z, carpeta)

    visualizacion(X, T, Z, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='plot_pequeño', cmap='seismic')
    plt.show()