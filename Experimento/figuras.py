from procesos import *
from visualizacion import *

if __name__ == '__main__':
    #[X, parametros, PHI_std, sinegauss_fit, gaussian_fit] = cargar_txt(experimental_data_file,
    #                                                                   data_file + '\datos_procesados', X='X',
    #                                                                   parametros='parametros', PHI_std='PHI_std',
    #                                                                   sinegauss_fit='sinegauss_fit',
    #                                                                   gaussian_fit='gaussian_fit')
    datos_path = 'F:\mnustes_science\experimental_data'
    carpeta = select_directory(datos_path)
    [X, T, Z] = cargar_txt(carpeta, '', X='X', T='T', Z='Z')
    visualizacion(X, T, Z, tipo='colormap', guardar='si', path=carpeta,
                  file='', nombre='plot_pequeño', cmap='seismic')