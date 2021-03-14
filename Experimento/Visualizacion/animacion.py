import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from IPython import display


if __name__ == '__main__':
    '''
    x = np.arange(-2, 2, 0.01)
    y = np.arange(-2, 2, 0.01)
    xx, yy = np.meshgrid(x, y)
    z = 1 * (xx ** 2 + yy ** 2)

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)

    data_skip = 10

    def init_func():
        ax.clear()
        plt.xlabel('pi')
        plt.ylabel('sin(pi)')

    def update_plot(i):
        print(i)
        ax.plot(xx[i:i + data_skip, :], z[i:i + data_skip, :], color='k')
        plt.xlim((x[0], x[-1]))
        plt.ylim((0, 1))

    anim = FuncAnimation(fig,
                         update_plot,
                         frames=np.arange(0, len(x), data_skip),
                         init_func=init_func,
                         interval=20)
    anim.save('basic_animation.gif', fps=30, writer='Pillow')
    '''

    x = np.linspace(0, 2 * np.pi, 100)

    fig = plt.figure()

    lines = plt.plot([])
    line = lines[0]
    #other setuo
    plt.xlim(0, 2*np.pi)
    plt.ylim(-1.1, 1.1)
    def animate(frame):
        y = np.sin(x + 2 * np.pi * frame/100)
        line.set_data((x, y))
        #update plot

    anim = FuncAnimation (fig, animate, frames=100, interval=20)
    anim.save('basic_animation.gif', fps=100)
    #video = anim.to_html5_video()
    #html = display.HTML(video)
    #display.display(html)
    plt.show()
