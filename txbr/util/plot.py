import numpy,pylab,matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial','Helvetica'],'size':'18'})
matplotlib.rc('lines',**{'lw':2,'markersize':6})
matplotlib.rc('axes',**{'titlesize':24})
matplotlib.rc('figure',**{'edgecolor':'w'})
matplotlib.rc('xtick',**{'labelsize':16})
matplotlib.rc('ytick',**{'labelsize':16})

def plot_contour(X,title='Contour'):
    '''Plot contours. Variable X is a 2D numpy array.
    '''

    n1,n2 = X.shape

    pylab.figure()
    pylab.contour(numpy.arange(n1),numpy.arange(n2),X.T)
    pylab.xlim(1,n1)
    pylab.ylim(1,n2)
    pylab.xlabel('X')
    pylab.ylabel('Y')
    pylab.title(title)