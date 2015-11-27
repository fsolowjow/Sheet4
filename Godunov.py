from __future__ import division
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from matplotlib import animation

def initial_value(x):
	if ( x < 0 ): return 0
	else:	return 1

# First set up the figure, the axis, and the plot element we want to animate
def set_up():
    global fig, ax, line, time_text
    fig = plt.figure()
    ax = plt.axes(ylim=(0, 1), xlim=(-2, 2))
    ax.set_ylabel('$q(x,t)$')
    ax.set_xlabel('$x$')
    ax.set_title('Godunovs method for linear systems')
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text,

# animation function.  This is called sequentially
def animate(j):
    y = U[:,j]
    time_text.set_text('Time = %.2f' % t[j])
    line.set_data(x,y)
    return line, time_text,
 
def run():
    set_up()
    # call the animator.  blit=True means only re-draw the parts that have change
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=range(t.size), interval=50, blit=False)
    plt.show()



a = 2

x_step = 0.05
t_step = 0.01

x = np.arange(-2, 2 , x_step)
t = np.arange ( 0 , 4 , t_step)

U = np.empty((x.size,t.size))

q = np.zeros( np.size(x) )
for i in range(np.size(q)) :
	q[i] = initial_value( x[i] )

for j in range(np.size(t)) :
	U[:,j] = q
	qtemp = q
	for i in range(np.size(q)):
		if( i == 0 ): qtemp[0] = qtemp[0] - ((a * t_step / x_step) * (qtemp[0] - qtemp[np.size(q) - 1]))
		qtemp[i] = qtemp[i] - ((a * t_step / x_step) * (qtemp[i] - qtemp[i-1]))
		if( i == np.size(q) - 1):	qtemp[i] = qtemp[i] - ((a * t_step / x_step) * (qtemp[i] - qtemp[0]))
	q = qtemp

run()