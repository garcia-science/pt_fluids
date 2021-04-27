
from procesos import *
from visualizacion import *

if __name__ == '__main__':
    ############ DETECTOR DE CONTORNOS ###########

    recs = ROI_select('E:\mnustes_science\images\img_lab\modo_4.1\\f=14.61_a=3.8')
    deteccion(images_file, '\img_lab' + '\modo_4.1\\f=14.61_a=3.8', '\canned' + '\modo_4.1\\f=14.61_a=3.8', recs, sigma=0.2)

    ######### GUARDADO Y CARGA DE DATOS #########
    IMGs = os.listdir(images_file + '\canned' + '\modo_4.1\\f=14.61_a=3.8')
    X, T, PHI = datos_3d(IMGs, images_file, '\canned' + '\modo_4.1\\f=14.61_a=3.8', nivel = 'si')
    guardar_txt(experimental_data_file, data_file, X=X, T=T, PHI=PHI)