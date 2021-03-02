import os
import numpy as np
from varname import nameof

path = r'C:\Users\mnustes_science\PT_fluids'
file_in = '\IMG'
file_out = '\IMG_canned'
reference_file = r'C:\Users\mnustes_science\PT_fluids\IMG'



def guardar_txt(array):
    for i in range(len(array)):
        np.savetxt("array0"+str(i)+".txt", array[i])

def cargar_txt(len):
    array = []
    for i in range(len):
        array_i = np.loadtxt("array0"+str(i)+".txt")
        array.append(array_i)
    return array