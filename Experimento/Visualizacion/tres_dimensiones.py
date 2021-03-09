from matplotlib import pyplot as plt
import numpy as np

def plot_ZXT(X, Y, Z):
    XX, YY = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig
def plot_XY(X, Y):
    plt.plot(X, Y)
    plt.xlabel(r'$x$ (Espacio)', fontsize=15)
    plt.ylabel(r'$\psi(x)$ (Altura)', fontsize=15)
    plt.title('Interesting Graph\nCheck it out')
    plt.grid(color='silver', linestyle='--', linewidth=0.5)

def color_map(X, Y, Z):
    XX, YY = np.meshgrid(X, Y)
    ax = plt.gca()
    pcm = ax.pcolormesh(X, Y, Z, vmin=np.min(Z), vmax=np.max(Z), cmap='RdBu_r', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label(r'        $\psi(x,t)$', rotation=0, fontsize=15)
    plt.xlabel(r'$x$ (Espacio)', fontsize=15)
    plt.ylabel(r'$t$ (Tiempo)', fontsize=15)
    plt.title('Interesting Graph\nCheck it out')
    return pcm

def DOS(frec, density, scale, titulo):
    plt.plot(frec, density)
    plt.yscale(str(scale))
    plt.xlabel(r'$\omega$ (Frecuencias)', fontsize=15)
    plt.ylabel(r'$\rho(\omega)$ (Densidad de Frecuencias)', fontsize=15)
    plt.title(str(titulo), fontsize=20)
    plt.grid(color='silver', linestyle='--', linewidth=0.5)