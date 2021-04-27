from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from directorios import *


def visualizacion(*args, tipo, guardar, path, file, nombre, **kwargs):
    if 'xlabel' not in kwargs:
        xlabel = ' '
    else:
        xlabel = kwargs['xlabel']
    if 'ylabel' not in kwargs:
        ylabel = ' '
    else:
        ylabel = kwargs['xlabel']
    if 'zlabel' not in kwargs:
        zlabel = ' '
    else:
        zlabel = kwargs['xlabel']
    if 'x_scale' not in kwargs:
        xscale = 'linear'
    else:
        xscale = kwargs['xscale']
    if 'yscale' not in kwargs:
        yscale = 'linear'
    else:
        yscale = kwargs['yscale']
    if 'titulo' not in kwargs:
        titulo = 'Agregar Título'
    else:
        titulo = kwargs['titulo']
    if 'color' not in kwargs:
        color = 'b'
    else:
        color = kwargs['color']
    if 'linestyle' not in kwargs:
        linestyle = '-'
    else:
        linestyle = kwargs['linestyle']
    if 'linewidth' not in kwargs:
        linewidth = 1
    else:
        linewidth = kwargs['linewidth']
    if 'cmap' not in kwargs:
        cmap = 'plasma'
    else:
        cmap = kwargs['cmap']
    if tipo == "2D":
        X = args[0]
        Y = args[1]
        plt.plot(X, Y, color=color, linestyle=linestyle, linewidth=linewidth)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        plt.xscale(xscale)
        plt.yscale(yscale)
        plt.title(titulo)
        plt.grid(color='silver', linestyle='--', linewidth=0.5)
        if guardar == 'si':
            plt.savefig(path + file + '\\' + nombre)
    elif tipo == "2D_multiple":
        X = args[0]
        Ys = args[1]
        if len(args) == 2:
            linestyles = ['-'] * len(Ys)
            colors = [(np.random.rand(3,))] * len(Ys)
        elif len(args) == 3:
            colors = [(np.random.rand(3,))] * len(Ys)
        else:
            colors = args[2]
            linestyles = args[3]
        for i in range(len(Ys)):
            plt.plot(X, Ys[i], linestyle=linestyles[i], color=colors[i])
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        plt.xscale(xscale)
        plt.yscale(yscale)
        plt.title(titulo)
        plt.grid(color='silver', linestyle='--', linewidth=0.5)
        if guardar == 'si':
            plt.savefig(path + file + '\\' + nombre)
    elif tipo == "colormap":
        X = args[0]
        Y = args[1]
        Z = args[2]
        ax = plt.gca()
        pcm = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), cmap=cmap, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label(zlabel, rotation=0, fontsize=15)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        plt.title(titulo)
        if guardar == 'si':
            plt.savefig(path + file + '\\' + nombre)
    elif tipo == "3D":
        X = args[0]
        Y = args[1]
        Z = args[2]
        XX, YY = np.meshgrid(X, Y)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False)
        plt.title(titulo)
        fig.colorbar(surf, shrink=0.5, aspect=5)
        if guardar == 'si':
            plt.savefig(path + file + '\\' + nombre)
    elif tipo == "animacion":
        campo = args[0]
        x = args[1]
        Nt = args[2]
        fig = plt.figure()
        lines = plt.plot([])
        line = lines[0]
        # other setuo
        plt.xlim(x[0], x[-1])
        plt.ylim(1.1 * np.amax(campo), 1.1 * np.amin(campo))
        plt.xlabel(kwargs['xlabel'], fontsize=15)
        plt.ylabel(kwargs['ylabel'], fontsize=15)
        plt.title(kwargs['titulo'])
        plt.grid(color='silver', linestyle='--', linewidth=0.5)
        def animate(frame):
            y = campo[100 * frame, :]
            line.set_data((x, y))
        anim = FuncAnimation(fig, animate, frames=int(Nt / 100), interval=60)
        if guardar == 'si':
            anim.save(path + file + '\\' + nombre + '.gif', fps=100)
        return anim