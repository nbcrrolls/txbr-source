import numpy,pylab
import util,util.plot
from util.plot import plot_contour

shape = (200,200)
points = numpy.array([[100,100,1.0],[100,150,-1.0],[100,300,1.0]])

map = util.mapSurface(shape,points)
util.plot.plot_contour(map.T,title='2D Interpolation test')
pylab.show()