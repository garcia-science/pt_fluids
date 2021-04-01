from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from Experimento.Input.enrutador import *


def func(x, a, b, c):
    return a * np.exp(-(x - b) ** 2 / (2 * c) ** 2)

def plot_ZXT(X, Y, Z, guardar, nombre, titulo):
    XX, YY = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    plt.title(titulo)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    if guardar == 'si':
        plt.savefig(plot_path + nombre)
    return fig


def plot_XY(X, Y, guardar, nombre, titulo, xlabel, ylabel):
    plt.plot(X, Y)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo)
    plt.grid(color='silver', linestyle='--', linewidth=0.5)
    if guardar == 'si':
        plt.savefig(plot_path + nombre)


def color_map(X, Y, Z, guardar, nombre, titulo, xlabel, ylabel):
    XX, YY = np.meshgrid(X, Y)
    ax = plt.gca()
    pcm = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), cmap='RdBu_r', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label(r'        $\psi(x,t)$', rotation=0, fontsize=15)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo)
    if guardar == 'si':
        plt.savefig(plot_path + nombre)
    return pcm


def DOS(frec, density, scale, guardar, nombre, titulo, xlabel, ylabel):
    plt.plot(frec, density)
    plt.yscale(str(scale))
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo, fontsize=20)
    plt.grid(color='silver', linestyle='--', linewidth=0.5)
    if guardar == 'si':
        plt.savefig(plot_path + nombre)


def animacion(campo, x_grid, Nt, L, y_lim, guardar, nombre, titulo, xlabel, ylabel):
    x = x_grid
    fig = plt.figure()
    lines = plt.plot([])
    line = lines[0]
    #other setuo
    plt.xlim(-L/2, L/2)
    plt.ylim(y_lim[0], y_lim[1])
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo)
    plt.grid(color='silver', linestyle='--', linewidth=0.5)
    def animate(frame):
        y = campo[100 * frame, :]
        line.set_data((x, y))
        #update plot
    anim = FuncAnimation(fig, animate, frames=int(Nt/100), interval=60)
    if guardar == 'si':
        anim.save(plot_path + nombre + '.gif', fps=100)
    return anim


def multiple_XY(X, Ys, guardar, nombre, titulo, xlabel, ylabel):
    for i in range(len(Ys)):
        plt.plot(X, Ys[i])
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo)
    plt.grid(color='silver', linestyle='--', linewidth=0.5)
    if guardar == 'si':
        plt.savefig(plot_path + nombre)


def color_map_for(X, Y, Z, guardar, path, file, nombre, titulo, xlabel, ylabel, z_name):
    XX, YY = np.meshgrid(X, Y)
    ax = plt.gca()
    pcm = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), cmap='RdBu_r', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label(z_name, rotation=0, fontsize=15)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(titulo)
    if guardar == 'si':
        plt.savefig(path + file + nombre)