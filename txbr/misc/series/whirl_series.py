import numpy, pylab, matplotlib
import modl

import matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial','Helvetica'],'size':'16'})
matplotlib.rc('lines',**{'lw':2,'markersize':6})
matplotlib.rc('axes',**{'titlesize':28})
matplotlib.rc('axes',**{'labelsize':22})
matplotlib.rc('figure',**{'edgecolor':'w'})
matplotlib.rc('xtick',**{'labelsize':16})
matplotlib.rc('ytick',**{'labelsize':16})

def plotXYZ(XYZ):

    fig = pylab.figure()

    ax = matplotlib.axes3d.Axes3D(fig)
    ax.scatter3D(XYZ[:,0],XYZ[:,1],XYZ[:,2])

    pylab.show()


def plotXY2(XYZ):

    fig = pylab.figure()

    pylab.scatter(XYZ[:,0],XYZ[:,1])


def plot_variations(r,radial_var,ortho_var,title):

    fig = pylab.figure()

    pylab.title(title)
    pylab.xlabel('Distance to Center')
    pylab.ylabel('Relative Variation')

    pylab.plot(r, radial_var, 'ob', label='Rel. radial Variation')
    pylab.plot(r, ortho_var, 'vr', label='Rel. orthoradial Variation')
    pylab.legend()

    pylab.ylim(-0.05,0.05)

    pylab.savefig(title.replace(" ","_") + '_var.png', format='png')


def plot_field(x,y,u,v,title):

    fig = pylab.figure()

    pylab.quiver(x, y, u, v, label='toto')

    pylab.title(title)
    pylab.xlabel('X')
    pylab.ylabel('Y')
    pylab.legend()

    pylab.savefig(title.replace(" ","_") + '_field.png', format='png')


def main(case,file):

    model = modl.Model(file)
    model.loadFromFile()

    object = model.getObject(0)

    points = object.points()

    size = points.shape[0]
    ntrack = object.numberOfContours()
    nz = size/ntrack

    print 'Number of Tracks=%i   Number of Steps=%i' %(ntrack,nz)

    # Fluctuations of tracks during the shooting

    d = numpy.zeros((ntrack),dtype=float)

    for iz in range(1,nz):
        d[:] = d[:] + (points[iz::nz,0]-points[0::nz,0])**2 \
            +  (points[iz::nz,1]-points[0::nz,1])**2 \
            +  (points[iz::nz,2]-points[0::nz,2])**2

    d = d/nz

    indexOfCenter = numpy.argmin(d)
    center = points[indexOfCenter*nz:(indexOfCenter+1)*nz,:]

    center = object.getContour(indexOfCenter).points

    print 'Index of Center; %i' %indexOfCenter
    print center

    xy = numpy.zeros((nz,ntrack-1,2),dtype=float)
    r = numpy.zeros((ntrack-1),dtype=float)
    radial_var = numpy.zeros((ntrack-1),dtype=float)
    ortho_var = numpy.zeros((ntrack-1),dtype=float)

    for itrack in range(1,ntrack):
        xy[:,itrack-1,0] = points[itrack*nz:(itrack+1)*nz,0]-center[:,0]
        xy[:,itrack-1,1] = points[itrack*nz:(itrack+1)*nz,1]-center[:,1]

    for itrack in range(ntrack-1):
        # define radial and orthoradial vector for step #0
        r_ = numpy.array([xy[0,itrack,0],xy[0,itrack,1]])
        u_ = numpy.array([xy[0,itrack,1],-xy[0,itrack,0]])
        r[itrack] = numpy.dot(r_,r_)**0.5
        radial_var[itrack] = numpy.dot(xy[-1,itrack,:]-xy[0,itrack,:],r_)/r[itrack]**2
        ortho_var[itrack] = numpy.dot(xy[-1,itrack,:]-xy[0,itrack,:],u_)/r[itrack]**2
        x_ = xy[0,itrack,0] + center[0,0]
        y_ = xy[0,itrack,1] + center[0,1]
        print 'Track #%i\t(%5.3f,%5.3f)\tr=%f\tradial=%5.3f\torthoradial=%f' \
            %(itrack,x_,y_,r[itrack],radial_var[itrack],ortho_var[itrack])

    x = xy[0,:,0]
    y = xy[0,:,1]
    u = (xy[-1,:,0]-xy[0,:,0])/r[itrack]
    v = (xy[-1,:,1]-xy[0,:,1])/r[itrack]

    plot_variations(r,radial_var,ortho_var,case)
    plot_field(x,y,u,v,case)

    #plotXY2(object.points())
    #plotXYZ(object.points())


if __name__ == '__main__':

    dict = { \
            'Stage Series 5k': '/Users/sph/Electron_Tomography/txbr/dataset/whirlpool/12_10_08_5k_20ctsStep_a.fid',\
            'Focus Series 5k':'/Users/sph/Electron_Tomography/txbr/dataset/whirlpool/12_10_08_5k_1micObjFoc.fid',\
            'Stage Series 10k':'/Users/sph/Electron_Tomography/txbr/dataset/whirlpool/12_10_08_10k_20ctsStep.fid',\
            'Focus Series 10k':'/Users/sph/Electron_Tomography/txbr/dataset/whirlpool/12_10_08_10k_1micObjFocReal.fid'\
            }

    for key in dict:
        main(key,dict[key])

    pylab.show()

