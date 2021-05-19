from procesos import *
from visualizacion import *
from matplotlib import pyplot as plt

if __name__ == '__main__':
    datos_path = 'D:\mnustes_science\experimental_data'
    root = tk.Tk()
    root.withdraw()
    carpeta = filedialog.askdirectory(parent=root,
                                      initialdir=datos_path,
                                      title='Selecciones la carpeta a procesar')
    [X, T, A, B, Z] = cargar_txt(carpeta, '', X='X', T='T', A='A', B='B', Z='Z')
    #plt.plot(T[0:200], Z[0:200, 670], '-o')
    #plt.plot(T[0:200], A[0:200, 670], 'r')
    #plt.plot(T[0:200], B[0:200, 670], 'g')
    #plt.grid(True)
    #plt.show()
    #A = filtro_superficie(A, 100, 'Y')
    #B = filtro_superficie(B, 100, 'Y')
    visualizacion(X, T, A, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='A_plot', cmap='seismic')
    plt.close()
    visualizacion(X, T, B, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='B_plot', cmap='seismic')
    plt.close()
