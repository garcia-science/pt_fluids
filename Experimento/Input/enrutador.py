import os
import numpy as np

path = r'C:\Users\mnustes_science\PT_fluids'
file_in = '\IMG'
file_out = '\IMG_canned'
reference_file = r'C:\Users\mnustes_science\PT_fluids\IMG'

def guardar_matriz(matriz, nombre):
    with open(nombre, 'wb') as f:
        for line in matriz:
            np.savetxt(f, line, fmt='%.2f')
def guardar_vector(vector, nombre):
    a_file = open(nombre, "w")
    np.savetxt(a_file, vector)
    a_file.close()