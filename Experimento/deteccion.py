
from procesos import *
from visualizacion import *


if __name__ == '__main__':
    ############ DETECTOR DE CONTORNOS ###########

    #os.makedirs(experimental_data_file + '\\07-05-2021\\f=15.22_a=6.0')
    recs = ROI_select('E:\mnustes_science\images\img_lab\\07-05-2021\\f=15.22_a=6.0')
    deteccion(images_file, '\img_lab' + '\\07-05-2021\\f=15.22_a=6.0', '\canned' + '\\07-05-2021\\f=15.22_a=6.0', recs, sigma=1)

    ######### GUARDADO Y CARGA DE DATOS #########

    IMGs = os.listdir(images_file + '\canned' + '\\07-05-2021\\f=15.22_a=6.0')
    X, T, PHI = datos_3d(IMGs, images_file, '\canned' + '\\07-05-2021\\f=15.22_a=6.0', nivel = 'si')
    guardar_txt(experimental_data_file, '\\07-05-2021\\f=15.22_a=6.0', X=X, T=T, PHI=PHI)
