import os
import numpy as np
from varname import nameof

path = r'C:\Users\mnustes_science\PT_fluids'
file_in = '\IMG'
file_out = '\IMG_canned'
reference_file = r'C:\Users\mnustes_science\PT_fluids\IMG'
data_path = r'C:\Users\mnustes_science\PT_fluids\pt_fluids\Experimento\datos'
data_file = r'\pruebas'
sim_path = r'C:\Users\mnustes_science\PT_fluids\pt_fluids\Simulaciones\pndls_sim'
plot_path = 'Plots'



def guardar_txt(array, path, file):
    for i in range(len(array)):
        np.savetxt(path + file + r"\datos0"+str(i)+".txt", array[i])


def cargar_txt(len, path, file):
    array = []
    for i in range(len):
        array_i = np.loadtxt(path + file + r"\datos0"+str(i)+".txt")
        array.append(array_i)
    return array


def nombre_pndls(gamma, mu, nu):
    gamma = gamma.split('.')
    mu = mu.split('.')
    nu = nu.split('.')
    nombre = '\gmn' + '_' + gamma[0] + gamma[1] + '_' + mu[0] + mu[1] + '_' + nu[0] + nu[1]
    return nombre