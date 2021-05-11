from tkinter import *
import tkinter as tk
from tkinter import filedialog
import os

from procesos import *
from visualizacion import *

if __name__ == '__main__':

    root = tk.Tk()
    root.withdraw()

    detection_parent_file = filedialog.askdirectory(parent=root,
                                                    initialdir="E:\mnustes_science",
                                                    title='Detección multiple')
    if not detection_parent_file:
        sys.exit('No se seleccionó ningún archivo')

    os.chdir(detection_parent_file)
    detection_files = os.listdir()
    parent_file_name = os.path.basename(detection_parent_file)
    canned_path = 'E:\mnustes_science\images\canned'
    datos_path = 'E:\mnustes_science\experimental_data'
    #ESTO ESTA ACA PARA TOMAR EL CERO DE REFERENCIA
    #reference_image = filedialog.askdirectory(parent=root,
    #                                                initialdir=detection_parent_file + '\\reference',
    #                                                title='Seleccionar imagen de referencia')

    for name in detection_files:
        reference_file = detection_parent_file + '\\' + name
        recs = ROI_select(reference_file)
        deteccion(detection_parent_file + '\\' + name, canned_path + '\\' + parent_file_name + '\\' + name, recs, sigma=1)
        IMGs = os.listdir(canned_path + '\\' + parent_file_name + '\\' + name)
        X, T, PHI = datos_3d(IMGs, canned_path + '\\' + parent_file_name + '\\' + name, nivel='si')
        guardar_txt(datos_path, '\\' + parent_file_name + '\\' + name, X=X, T=T, PHI=PHI)
