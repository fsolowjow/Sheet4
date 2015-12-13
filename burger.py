from __future__ import division
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from matplotlib import animation

#-------------------------------------------Begin animation stuff----------------------------------------------	
# First set up the figure, the axis, and the plot element we want to animate
def set_up():
    global fig, ax, line, time_text
    fig = plt.figure()
    ax = plt.axes(ylim=(-2, 2), xlim=(-2, 2))
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
#-------------------------------------------End animation stuff-------------------------------------------------

	
def initial_cond1(x):
	return np.exp(-4 * (x-1)*(x-1))

def initial_cond2(x):
	if( isinstance(x, float) == True or isinstance(x, int) == True):
		if( x > 0 ):	return 0
		else:			return 1
	else:
		temp = np.zeros(np.size(x))
		for i in range(np.size(x)):
			if( x[i] <= 0 ):	temp[i]=1
		return temp
	
def initial_cond3(x):
	if( isinstance(x, float) == True or isinstance(x, int) == True):
		if( x <= 0):			return 0
		if( x < 1 and x > 0):	return (1-x)
		if(	x >= 1):			return 1
	else:
		temp = np.zeros(np.size(x))
		for i in range(np.size(x)):
			if( x[i] < 1 and x[i] > 0):	temp[i] = (1-x[i])
			if(	x[i] >= 1):				temp[i] = 1
		return temp
		
def initial_cond4(x):
	if( isinstance(x, float) == True or isinstance(x, int) == True):
		if( x > 0 ):	return 1
		else:			return 0.5
	else:
		temp = np.zeros(np.size(x))
		for i in range(np.size(x)):
			if( x[i] <= 0 ):	temp[i]=0.5
			else:				temp[i]=1
		return temp
		
def initial_cond (x,i):
	if(i==1): return initial_cond1(x)
	if(i==2): return initial_cond2(x)
	if(i==3): return initial_cond3(x)	
	if(i==4): return initial_cond4(x)	

def update_q ( q , x_step , t_step):
	qtemp = q
	for i in range( np.size(q)):
		if ( i == 0 ):							q[i] = qtemp[i] - (t_step / x_step) * ( flow(qtemp[i] , qtemp[i+1]) - flow(qtemp[np.size(q)-1],qtemp[i]) )
		if ( i == np.size(q)-1) :				q[i] = qtemp[i] - (t_step / x_step) * ( flow(qtemp[i] , qtemp[0]) - flow(qtemp[i-1],qtemp[i]) )
		if ( i!= 0 and i != np.size(q)-1 ):		q[i] = qtemp[i] - (t_step / x_step) * ( flow(qtemp[i] , qtemp[i+1]) - flow(qtemp[i-1],qtemp[i]) )
	return qtemp
		
def flow ( q1 , q2 ):
	temp=0
	if(q1 < q2 ):
		if ( q1 > 0):				temp = q1
		if ( q2 < 0):				temp = q2
		if ( q1 <= 0 and q2 >= 0):	temp = 0
	else:
		if( (q1+q2)/2 > 0 ):		temp = q1
		else:						temp = q2
	return temp*temp/2
	
def burgers_solver(x_l, x_r, t_end , eta, initial_condition , x_step):
						
	q = np.zeros( np.size(x))
	q = initial_cond( x , initial_condition)
	
	U = q
	t_step = eta * x_step / max(abs(q))
	t_counter = 0
	
	while( t_counter < t_end):
		q = update_q(q , x_step , t_step)
		t_counter += t_step
		U = np.vstack((U , q))

	return U

x_l = input(" Value for x_l: (-2)\n")
x_r = input(" Value for x_r: (2)\n")
t_end = input(" Value for t_end: (4)\n")
eta = input(" Value for eta: (0.9)\n")

x_step = 0.01
x = np.arange(x_l, x_r , x_step)

initial_condition = input(" Initial conditions: \n 1 Gaussian \n 2 Riemann \n 3 Kink \n 4 Other Riemann ")
		
U = burgers_solver(x_l, x_r, t_end , eta, initial_condition , x_step)
#saved values in U the wrong way around, because vstack is strange => U.T	
U= U.T

t = np.arange ( 0 , t_end , t_end / np.size(U[0,:]))

run()
