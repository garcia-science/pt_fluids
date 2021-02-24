from matplotlib import pyplot as plt
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *

if __name__ == '__main__':
    #recs = ROI_select(reference_file)
    #deteccion(path, file_in, file_out, recs, 0.37)
    IMGs = os.listdir(path + file_out)
    PHI, T, X, Y, Z = datos_3d(IMGs, path, file_out, filtro = 'si')

    #guardar_matriz(PHI, nombre = "PHI.txt")
    #guardar_vector(X, nombre = "X.txt")

    #PHI = np.loadtxt('PHI.txt')
    #X = np.loadtxt('X.txt')
    #fig = plot_ZXT(X, Y, Z)
    #PHIT_proy = proyeccion(Z)
    #fig1 = plt.plot(X, PHIT_proy)
    fig2 = color_map(X, Y, Z)
    #plot = plt.plot(X, Z[0,:])
    plt.show()