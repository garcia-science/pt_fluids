from matplotlib import pyplot as plt
import numpy as np

def plot_ZXT(XX, YY, Z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig