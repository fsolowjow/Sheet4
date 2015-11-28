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
    ax = plt.axes(ylim=(0, 2), xlim=(-2, 2))
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


#Decomposes initial data into the components
def initial_values(q_l , q_r , eigenvector , x):
	dim = np.size(q_l)
	qtemp = np.zeros( (np.size(x) , dim) )
	w_r = LA.solve(eigenvector, q_r)
	w_l = LA.solve(eigenvector, q_l)
	for i in range(dim):
		for j in range(np.size(x)):
			if(x[j] < 0):	qtemp[j,i] = w_l[i]
			else:			qtemp[j,i] = w_r[i]
	
	return qtemp
	
def initial_values_1dim(x):
	if ( x < 0 ): return 0
	else:	return 1

def Godunov_linear_solv(A,q_l,q_r):

	dim = np.size(q_l)
	#Distinguish between 1dim case and system.
	#Case of a system: 
	if(dim > 1):
		eigenvalue , eigenvector = LA.eig(A)
		U = np.empty((x.size,t.size))
		#Sets the Q_i up for the first time step with the initial data. 
		#Q[j,i] is a matrix and contains the values for component i at x[j] at each time step
		q = initial_values( q_l , q_r , eigenvector , x )

		#iterating over time
		for j in range(np.size(t)) :
			#Values for the animation are saved in U
			U[:,j] = q[:,2]									#Change here to animate other components
			#The values for the next time step are saved in qtemp and afterwards q -> qtemp
			qtemp = q
			for n in range( dim ):
				for i in range(np.size(x)):
					if( i == 0 ): qtemp[0,n] = qtemp[0,n] - (( eigenvalue[n] * t_step / x_step) * (qtemp[0,n] - qtemp[np.size(x) - 1 ,n]))
					qtemp[i,n] = qtemp[i,n] - ((eigenvalue[n] * t_step / x_step) * (qtemp[i,n] - qtemp[i-1,n]))
					if( i == np.size(x) - 1):	qtemp[i,n] = qtemp[i,n] - ((eigenvalue[n] * t_step / x_step) * (qtemp[i,n] - qtemp[0,n]))
				q = qtemp
	#1dim case with the wavespeed A, that gets passed instead of a matrix:			
	else:
		U = np.empty((x.size,t.size))
		q = np.zeros( np.size(x) )
		for i in range(np.size(q)) :
			q[i] = initial_values_1dim( x[i] )

		for j in range(np.size(t)) :
			U[:,j] = q
			qtemp = q
			for i in range(np.size(q)):
				if( i == 0 ): qtemp[0] = qtemp[0] - ((A * t_step / x_step) * (qtemp[0] - qtemp[np.size(q) - 1]))
				qtemp[i] = qtemp[i] - ((A * t_step / x_step) * (qtemp[i] - qtemp[i-1]))
				if( i == np.size(q) - 1):	qtemp[i] = qtemp[i] - ((A * t_step / x_step) * (qtemp[i] - qtemp[0]))
			q = qtemp
	
	return U

t_step = 0.01
#t_step = input("Enter time stepsize:")
t = np.arange ( 0 , 4 , t_step)

x_step = 0.05
#x_step = input("Enter space stepsize:")
x = np.arange(-2, 2 , x_step)

A3 = np.array ( [[2 , 1 , -1] , [1, 1 , 1] , [-1, 1 , 2]])
q_l3 = np.array([1,0,1])
q_r3 = np.array([0,1,1])	

A2 = np.array ( [[2 , 1] , [0.0001,2]])
q_l2 = np.array([0,1])
q_r2 = np.array([1,0])	

A1 = 2
#q not realy needed for 1dim case
q_l1 = 0
q_r1 = 1

U=Godunov_linear_solv( A3 , q_l3 , q_r3)
run()

# The programm can crash for max eigenvalue * c > 1, because the CFL condition is not satisfied and therefore it gets instable