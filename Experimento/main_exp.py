from matplotlib import pyplot as plt

from Experimento.Visualizacion.tres_dimensiones import plot_ZXT
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *

if __name__ == '__main__':
    #recs = ROI_select(reference_file)
    #deteccion(path, file_in, file_out, recs, 4)

    PHI, T, XX, YY, Z = datos_3d(IMGs, path, file_out, 'si')
    fig = plot_ZXT(XX, YY, Z)
    plt.show()