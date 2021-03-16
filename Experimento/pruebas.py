from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
from Experimento.Visualizacion.tres_dimensiones import *
from Input.enrutador import *
from Source.deteccion import *
from Source.procesos import *
import cv2
import os
from scipy.optimize import curve_fit

if __name__ == '__main__':
    def func(x, a, b, c):
        return a * np.exp(-(x - b) ** 2 / (2 * c) ** 2)
    xdata = np.linspace(0, 4, 50)
    y = func(xdata, 2.5, 1.3, 0.5)
    np.random.seed(1729)
    y_noise = 0.2 * np.random.normal(size=xdata.size)
    ydata = y + y_noise
    plt.plot(xdata, ydata, 'b-', label='data')
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata, func(xdata, *popt), 'r-',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()