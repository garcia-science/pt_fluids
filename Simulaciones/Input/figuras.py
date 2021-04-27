from procesos import *
from visualizacion import *

if __name__ == '__main__':
    [X, parametros, PHI_std, sinegauss_fit, gaussian_fit] = cargar_txt(experimental_data_file,
                                                                       data_file + '\datos_procesados', X='X',
                                                                       parametros='parametros', PHI_std='PHI_std',
                                                                       sinegauss_fit='sinegauss_fit',
                                                                       gaussian_fit='gaussian_fit')
    visualizacion(X, [PHI_std, sinegauss_fit, gaussian_fit], ['r', 'grey', 'black'], ['-', ':', '--'], tipo="2D_multiple", guardar='si',
                  path=experimental_data_file, file=data_file + '\datos_procesados',
                  nombre='plot')