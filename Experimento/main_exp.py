from matplotlib import pyplot as plt
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *

if __name__ == '__main__':
    #recs = ROI_select(reference_file)
    #deteccion(path, file_in, file_out, recs, 3)
    IMGs = os.listdir(path + file_out)
    PHI, T, X, Y, Z = datos_3d(IMGs, path, file_out, 'si')
    #fig = plot_ZXT(X, Y, Z)
    fig = color_map(X, Y, Z)
    plt.show()