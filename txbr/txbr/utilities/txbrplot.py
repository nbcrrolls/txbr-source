import os
import sys
import numpy
import matplotlib
import pylab

matplotlib.rc('text', usetex=True)
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial','Helvetica'],'size':'18'})
matplotlib.rc('lines',**{'lw':2,'markersize':6})
matplotlib.rc('axes',**{'titlesize':24})
matplotlib.rc('figure',**{'edgecolor':'w'})
matplotlib.rc('xtick',**{'labelsize':16})
matplotlib.rc('ytick',**{'labelsize':16})


def show():

    pylab.show()
    

def plot3DMarkers( X, Y, Z, scale=(1.0,1.0,1.0), surfaces=() ):
    '''Surfaces is a list containing points of a surface'''
    
    try:
        from enthought.mayavi import mlab
    except ImportError:
        log.info('Mlab is not available!')
        return

    scale = numpy.array(scale)/25.0
    
    x = X.astype(float).ravel()/scale[0]
    y = Y.astype(float).ravel()/scale[1]
    z = Z.astype(float).ravel()/scale[2]
    
    fig = mlab.figure()
    
    mlab.points3d(x,y,z,scale_factor=1.0)
    
    for s in surfaces:
        x = s[0].astype(float)/scale[0]
        y = s[1].astype(float)/scale[1]
        z = s[2].astype(float)/scale[2]
        mlab.surf(x,y,z)
        
    mlab.show()
        
        
def saveAndLinkPlot( directory, basename, extension, model ):
    
    if directory==None or basename==None or not os.path.exists(directory):
        return

    try:
    
        current_file = '%s_%s_%s.png' %( basename, model, extension )
        current_file = os.path.join( directory, current_file )

        pylab.savefig( current_file, format='png' )

        file = '%s_%s.png' %( basename, extension )
        file = os.path.join( directory, file )

        if os.path.lexists(file): os.unlink(file)

        os.symlink(os.path.basename(current_file),file)
    
    except:

        print "Unexpected error in saving plot:", sys.exc_info()[0]
         
        
def plotErrorByTilt( E1byTilt, E2byTilt, EbyTilt, directory=None, basename=None, model=""):

    ntilt = EbyTilt.size
    tilts = numpy.arange(1,ntilt+1)

    pylab.figure()

    pylab.plot(tilts,EbyTilt,'bo-',tilts,E1byTilt, 'r:d',tilts,E2byTilt, 'g:1')

    pylab.xlim(1,ntilt)

    pylab.legend(('E',r'$E_P$',r'$E_T$'))
    pylab.title(r'Square Distance Error by Tilt')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'Error')
    
    saveAndLinkPlot( directory, basename, "err-by-tilt", model )

    
def plotErrorByPatch( E1byPatch, E2byPatch, EbyPatch, directory=None, basename=None, model=""):

    npatch = EbyPatch.size
    patches = numpy.arange(1,npatch+1)

    pylab.figure()

    pylab.plot(patches,EbyPatch,'bo-',patches,E1byPatch, 'r:d',patches,E2byPatch, 'g:1')

    pylab.xlim(1,npatch)

    pylab.legend(('E',r'$E_P$',r'$E_T$'))
    pylab.title(r'Square Distance Error by Patch')
    pylab.xlabel(r'Patch \#')
    pylab.ylabel(r'Error')
    
    saveAndLinkPlot( directory, basename, "err-by-patch", model )
    
    
def plotErrorByTiltAndPatch( res1, res2, res, directory=None, basename=None, model=""):

    ntilt,npatch = res.shape
    tilts = numpy.arange(1,ntilt+1)

    pylab.figure()

    pylab.plot(tilts,res,'bo-',tilts,res1, 'r:d',tilts,res2, 'g:1')

    pylab.xlim(1,ntilt)

    pylab.legend(('e',r'$e_{P}$',r'$e_{T}$'))
    pylab.title(r'Square Distance Error by Tilt and Patch')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'Error')
    
    saveAndLinkPlot( directory, basename, "err-by-tilt-and-patch", model )
    

def plot_projections( P, P_id, directory=None, basename=None, model="", save=True, show=False ):
    '''This routine plots the projection coefficients as a function of the tilt
    index'''
    
    (dim,nterm,ntilt) = P.shape

    t = range(ntilt)

    # Translation Coefficients

    pylab.figure()
    
    pylab.plot(P[0,0,:],'b')
    pylab.plot(P[1,0,:],'r')
    pylab.scatter(t,P_id[0,0,:],c='b',marker='o',s=15)
    pylab.scatter(t,P_id[1,0,:],c='r',marker='h',s=15)

    pylab.xlim(0,ntilt)
    pylab.legend((r'$t_{1}$',r'$t_{2}$'))
    pylab.title(r'Projection Coefficients')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'Translation Coefficients')
    
    if save: saveAndLinkPlot( directory, basename, "tr-coeffs", model )

    # X coefficients

    pylab.figure()
    pylab.plot(t,P[0,1,:],'b')
    pylab.plot(t,P[0,2,:],'r')
    pylab.plot(t,P[0,3,:],'g')
    pylab.scatter(t,P_id[0,1,:],c='b',marker='o',s=15)
    pylab.scatter(t,P_id[0,2,:],c='r',marker='h',s=15)
    pylab.scatter(t,P_id[0,3,:],c='g',marker='^',s=15)

    pylab.xlim(0,ntilt)
    pylab.legend((r'$g_{11}$',r'$g_{12}$',r'$g_{13}$'))
    pylab.title(r'Projection Coefficients')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'x linear coefficients')
    
    if save: saveAndLinkPlot( directory, basename, "x-coeffs", model )

    # Y coefficients

    pylab.figure()
    pylab.plot(t,P[1,1,:],'b')
    pylab.plot(t,P[1,2,:],'r')
    pylab.plot(t,P[1,3,:],'g')
    pylab.scatter(t,P_id[1,1,:],c='b',marker='o',s=15)
    pylab.scatter(t,P_id[1,2,:],c='r',marker='h',s=15)
    pylab.scatter(t,P_id[1,3,:],c='g',marker='^',s=15)

    pylab.xlim(0,ntilt)
    pylab.legend((r'$g_{21}$',r'$g_{22}$',r'$g_{23}$'))
    pylab.title(r'Projection Coefficients')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'y linear coefficients')

    if save: saveAndLinkPlot( directory, basename, "y-coeffs", model )

    if show: pylab.show()


def plot_reprojection( xy, xy_r, mask, directory=None, basename=None, model="", limits=None ):
    '''This routine is designed to plot the fiducial traces and their
    reprojection.
    Variable xy: a numpy array of shape (ntilt,npatch,n,2) that represents the 
    original marker traces.
    Variable xy_r: a numpy array of shape (ntilt,npatch,n,2) that represents the 
    marker reprojection
    Variable mask: indicates if a data point should be used for a given tilt and
    a given marker
    '''
    
    ntilt = mask.shape[0]
    npatch = mask.shape[1]

    ntilt_eff = numpy.sum(numpy.any(mask,axis=1))
    npatch_eff = numpy.sum(numpy.any(mask,axis=0))
    
#    print "Effective Number of tilt in the plot: %i" %ntilt_eff
#    print "Effective Number of patch int the plot: %i" %npatch_eff
    
    indexOfTilt = None
    indexOfPatch = None
    
    if ntilt_eff==1: 
        indexOfTilt = numpy.where(numpy.any(mask,axis=1))[0]
        
    if npatch_eff==1: 
        indexOfPatch = numpy.where(numpy.any(mask,axis=0))[0]

    pylab.figure()
    
    for ipatch in range(npatch):
        indices = numpy.where(mask[:,ipatch])
        if len(indices[0])==0: 
            continue
        xy_ = xy[:,ipatch,:][indices]
        xy_r_ = xy_r[:,ipatch,:][indices]
        pylab.scatter( xy_[:,:,0], xy_[:,:,1], c = 'green', label='Tracks' )
        pylab.plot( xy_r_[:,:,0], xy_r_[:,:,1], 'r:1', label='Reprojected Data')

    pylab.legend(('Reprojected Data',))

    #pylab.title(r'Contour Reprojection')
    pylab.xlabel(r'x')
    pylab.ylabel(r'y')

    if limits!=None:
        pylab.axis(limits)

    if indexOfTilt!=None and indexOfPatch!=None:
        extension = "_tilt-%i_patch-%i" %(indexOfTilt+1,indexOfPatch+1)
        pylab.figtext(0.15,0.85,r'Patch \#%i/%i  Tilt \#%i/%i' %(indexOfPatch+1,npatch,indexOfTilt+1,ntilt),fontsize=14)
    elif indexOfTilt!=None:
        extension = "_tilt-%i" %(indexOfTilt+1)
        pylab.figtext(0.15,0.85,r'Tilt \#%i/%i' %(indexOfTilt+1,ntilt),fontsize=14)
    elif indexOfPatch!=None:
        extension = "_patch-%i" %(indexOfPatch+1)
        pylab.figtext(0.15,0.85,r'Patch \#%i/%i' %(indexOfPatch+1,npatch),fontsize=14)
    else:
        extension = ""
        
    saveAndLinkPlot( directory, basename, "reproj", model )


def plot_angles( angles, anglesid ):
    '''This routine is used to plot the angles intervening in the orthogonal 
    approximation for the TxBR projection map. They are plotted versus the 
    tilt index.'''
    
    (dim,ntilt) = angles.shape
    
    t = range(ntilt)
    
    pylab.figure()
    
    pylab.plot(t,angles[0,:],'b')
    pylab.plot(t,angles[1,:],'r')
    pylab.plot(t,angles[2,:],'g')
    pylab.scatter(t,anglesid[0,:],c='b',marker='o',s=8)
    pylab.scatter(t,anglesid[1,:],c='r',marker='h',s=8)
    pylab.scatter(t,anglesid[2,:],c='g',marker='^',s=8)

    pylab.xlim(0,ntilt)
    pylab.legend((r'$\alpha$',r'$\beta$',r'$\gamma$'))
    pylab.title(r'Orthogonal Coefficients')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'$\alpha, \beta, \gamma$')
    
    
def plot_magnifications( magnification ):
    '''This routine is used to plot the magnification factors intervening in the orthogonal 
    approximation for the TxBR projection map. They are plotted versus the 
    tilt index.'''
    
    ( dim, ntilt ) = magnification.shape
    
    t = range(ntilt)
    
    pylab.figure()
    
    for index in range(dim):
        pylab.plot(t,magnification[index,:])
        
    pylab.scatter(t,numpy.ones((ntilt)),marker='o',s=8)

    pylab.xlim(0,ntilt)
    pylab.title(r'Magnification Factor(s)')
    pylab.xlabel(r'Tilt \#')
    pylab.ylabel(r'M')
    
    
def plot_contours( pos1, delta, inputs, outputs, nx, ny, npoints):

    pos1 = numpy.resize(pos1,(pos1.size/npoints,npoints))
    delta = numpy.resize(delta,(delta.size/npoints,npoints))

    pos2 = pos1 - delta

    pos1 = numpy.swapaxes(pos1,0,1)
    pos2 = numpy.swapaxes(pos2,0,1)

    inp = numpy.row_stack(inputs)
    out = numpy.row_stack(outputs)

    pylab.figure()

    pylab.plot(pos1.real,pos1.imag,'b')
    pylab.plot(pos2.real,pos2.imag,'r')

    pylab.scatter(inp[:,0],inp[:,1],c='b')
    pylab.scatter(out[:,0],out[:,1],c='r')

    pylab.xlim(0,nx)
    pylab.ylim(0,ny)
    pylab.xlabel('X')
    pylab.ylabel('Y')
    
    
def plotRemap( data_in, data_out, directory=None, basename=None, model="" ):
    
    ntilts = len(data_in)
    
    pylab.figure()
    
    for itilt in range(ntilts):
        
        data_in_ = data_in[itilt]
        data_out_ = data_out[itilt]
        
        if data_in_.shape[0]!=0:
            pylab.scatter(data_in_[:,0],data_in_[:,1],c='b',marker='x')
        if data_out_.shape[0]!=0:
            pylab.scatter(data_out_[:,0],data_out_[:,1],c='r')
        
    pylab.xlabel("X")
    pylab.ylabel("Y")
    pylab.title("2D remap")
    
    saveAndLinkPlot( directory, basename, "remap", model )
    
    
def plotGaugeCharacterizarion( ellipses, tgts_x, limits, directory=None, basename=None, model="" ):
    
    xmin,xmax,ymin,ymax,zmin,zmax = limits
    
    pylab.figure()
        
    plotTangentLines(tgts_x,[ymin,ymax])
        
    for ellipse_left,ellipse_right in ellipses:
        plotEllipse(ellipse_left,'left')
        plotEllipse(ellipse_right,'right')
            
    pylab.xlim(xmin,xmax)
    pylab.ylim(ymin,ymax)
    
    if basename!=None:
        basename = basename.replace("_","\_")
        pylab.title(basename)
    
    saveAndLinkPlot( directory, basename, "ellipses", model )
     
    
def plotTangentLines(tangents,zlim):
    ''' Print the envelope.'''

    zlim = numpy.asarray(zlim)

    one = numpy.ones((2))
    t = numpy.row_stack((one,zlim))

    s = tangents.shape

    xpts = - numpy.tensordot(tangents[:,:,[0,2]],t,([2],[0]))

    alpha = tangents[:,:,1].repeat(2)
    alpha = numpy.resize(alpha,(s[0],s[1],2))
    xpts = xpts/alpha

    size = s[0]*s[1]

    xpts = numpy.resize(xpts,(size,2))
    zpts = numpy.tile(zlim,(size,1))

    pylab.plot(xpts.T,zpts.T,'r:')
        
    
def plotEllipse( ellipse, option=None ):
    '''Plot an ellipse'''
    
    x0 = ellipse.cx
    y0 = ellipse.cy
    a = ellipse.a
    b = ellipse.b
    theta = ellipse.phi

    if abs(a)<abs(b):
        tmp = a
        a = b
        b = tmp
        theta = theta - numpy.sign(theta)*numpy.pi/2

    theta0 = 5*numpy.pi/12;

    theta_1 = 0;
    theta_2 = 2*numpy.pi;

    epsilon = numpy.sqrt(1-b**2/a**2)

    if option=='left':
        r0 = a*(1-epsilon**2)/(1 + epsilon*numpy.cos(theta0))
        theta_1 = numpy.pi - numpy.arctan(r0*numpy.sin(theta0)/(r0*numpy.cos(theta0)+2*epsilon*a)) - theta
        theta_2 = numpy.pi + numpy.arctan(r0*numpy.sin(theta0)/(r0*numpy.cos(theta0)+2*epsilon*a)) - theta

    if option=='right':
        theta_1 = - theta0 - theta
        theta_2 = theta0 - theta

    focal_x = x0 + epsilon*a*numpy.cos(theta);
    focal_y = y0 + epsilon*a*numpy.sin(theta);

    phi = numpy.linspace(theta_1,theta_2,num=50)

    r = a*(1-epsilon**2)/(1+epsilon*numpy.cos(phi))

    X = focal_x + r*numpy.cos(phi+theta)
    Y = focal_y + r*numpy.sin(phi+theta)

    pylab.xlabel('X')
    pylab.ylabel('Z')
    pylab.plot(X,Y,'b-')
    

def plot_deviation_field( map0, map, stepx=25, stepy=25):
    '''This routine plots a 2D deviation field.
    map0: complex array that represents an undeformed grid
    map: complex array that represents an deformed grid
    '''

    nx,ny = map0.shape

    u = map0 - map

    pylab.figure()
    pylab.contour(numpy.arange(ny),numpy.arange(nx),u.real)

    pylab.figure()
    pylab.contour(numpy.arange(ny),numpy.arange(nx),u.imag)

    x0 = map0.real[::stepx,::stepy]
    y0 = map0.imag[::stepx,::stepy]

    ux = u.real[::stepx,::stepy]
    uy = u.imag[::stepx,::stepy]

    u = numpy.sqrt(ux**2+uy**2)

    ux = ux/numpy.sqrt(u)
    uy = uy/numpy.sqrt(u)

    pylab.figure()
    pylab.quiver(x0,y0,ux,uy,scale=150.0)
    pylab.xlabel('X')
    pylab.ylabel('Y')
    pylab.title('Displacement Field')


def plot_map_file(filename):

    if not os.path.exists(filename):
        print 'File %s does not exist!' %(filename)
        return

    f = mrc.MRCFile(filename)
    map = f.getZSliceAt(0)

    x = numpy.arange(f.nx,dtype='float32')
    y = numpy.arange(f.ny,dtype='float32')

    map0 = transpose(*pylab.meshgrid(x,y))

    plot_deviation_field(map0,map)
    