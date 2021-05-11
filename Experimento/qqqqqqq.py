import tkinter as tk
from tkinter import filedialog
import os

from procesos import *
from visualizacion import *

if __name__ == '__main__':
    print('hola')
    root = tk.Tk()
    root.withdraw()
    detection_parent_file = filedialog.askopenfilename(parent=root,
                                                    initialdir="E:\mnustes_science",
                                                    title='Detección multiple')
    print(detection_parent_file)