import os.path,time,math
import logging,logging.config
import numpy,matplotlib
import scipy,scipy.interpolate,pylab
import cv
import mrc,modl,util
import txbr.utilities

ISO_VOLUME_WIDTH = 'isov'
MINIMUM_WIDTH = 'min'
MAXIMUM_WIDTH = 'max'

POL_APPR_MODE = 'polynom'
SPLINE_APPR_MODE = 'spline'
INTERP_APPR_MODE = 'interpolation'

X_INTERP = 256    # Interpolation method is always done for the same X_width
Y_INTERP = 256    # Interpolation method is always done for the same Y_width

log = logging.getLogger()

def permut(file,direction):

    if direction=='X':
        n1 = file.ny
        n2 = file.nx
        n3 = file.nz
        slicaxis = 'Z'
    elif direction=='Y':
        n1 = file.nz
        n2 = file.ny
        n3 = file.nx
        slicaxis = 'X'
    else:
        n1 = file.nx
        n2 = file.nz
        n3 = file.ny
        slicaxis = 'Y'

    return (n1,n2,n3,slicaxis)

def eval_spl_boundaries(trk_path,direction='Z',mode=POL_APPR_MODE,order=[10,10],doPlot=True):
    ''' This routine approximates the external specimen surfaces with either polynomial functions
    (mode='polynom') , splines (mode='spline'), or with another interpolation-like approach. Variable
    'direction' indicates the axis variable that is function of the two others; it roughly corresponds
    to the specimen normal. Variable 'order' is a 2D numpy-array containing the maximal order in the
    free variables used for the polynomial regression.
    Returns (box,spec_limits,approximation)
    - box: the volume boundaries
    - spec_limits: the minimal box containing the specimen (or the points used for the boundaries)
    - approximation: the calculated approximation elements
    '''

    if trk_path==None or not os.path.exists(trk_path):
        raise IOError, 'File %s does not exists!' %trk_path

    model = modl.Model(trk_path)
    model.loadFromFile()

    # Tracked points on both surfaces

    pts_1_ = model.objects[0].points()
    pts_2_ = model.objects[1].points()

    # Keep going

    pts_1 = numpy.zeros(pts_1_.shape)
    pts_2 = numpy.zeros(pts_2_.shape)

    if direction=='X':
        xmax = model.ymax
        ymax = model.zmax
        zmax = model.xmax
        pts_1[:,0] = pts_1_[:,1]
        pts_1[:,1] = pts_1_[:,2]
        pts_1[:,2] = pts_1_[:,0]
        pts_2[:,0] = pts_2_[:,1]
        pts_2[:,1] = pts_2_[:,2]
        pts_2[:,2] = pts_2_[:,0]
    elif direction=='Y':
        xmax = model.zmax
        ymax = model.xmax
        zmax = model.ymax
        pts_1[:,0] = pts_1_[:,2]
        pts_1[:,1] = pts_1_[:,0]
        pts_1[:,2] = pts_1_[:,1]
        pts_2[:,0] = pts_2_[:,2]
        pts_2[:,1] = pts_2_[:,0]
        pts_2[:,2] = pts_2_[:,1]
    else:
        xmax = model.xmax
        ymax = model.ymax
        zmax = model.zmax
        pts_1 = pts_1_
        pts_2 = pts_2_

    box = [[1,xmax],[1,ymax],[1,zmax]]

    if (min(pts_1[:,2])>min(pts_2[:,2])):
        tmp = pts_1
        pts_1 = pts_2
        pts_2 = tmp

    min_1 = numpy.min(pts_1,axis=0)
    max_1 = numpy.max(pts_1,axis=0)
    min_2 = numpy.min(pts_2,axis=0)
    max_2 = numpy.max(pts_2,axis=0)

    range_init_1 = max_1[2]-min_1[2]
    range_init_2 = max_2[2]-min_2[2]

    spec_limits = \
            [[min(min_1[0],min_2[0]),max(max_1[0],max_2[0])],
             [min(min_1[1],min_2[1]),max(max_1[1],max_2[1])],
             [min(min_1[2],min_2[2]),max(max_1[2],max_2[2])]]

    log.info('Initial Extension of the boundary #1: %.2f' %range_init_1)
    log.info('                                  #2: %.2f' %range_init_2)

    if mode==POL_APPR_MODE:    # Polynomial Regression of the boundary surfaces

        order = numpy.asarray(order)

        log.info('Polynomial Regression with order %s!' %(order))

        parameters_1,base_1 = pts_1[:,:2],pts_1[:,2]
        poly_reg_1 = util.polynomial_regression(parameters_1,base_1,order)

        parameters_2,base_2 = pts_2[:,:2],pts_2[:,2]
        poly_reg_2 = util.polynomial_regression(parameters_2,base_2,order)

        poly_reg = [poly_reg_1,poly_reg_2]

        diff1 = poly_reg_1.eval(parameters_1)-base_1
        diff1_range = max(diff1)-min(diff1)

        diff2 = poly_reg_2.eval(parameters_2)-base_2
        diff2_range = max(diff2)-min(diff2)

        log.info('Boundary #1: %s' %(poly_reg_1))
        log.info('Boundary #2: %s' %(poly_reg_2))
        log.info('Deviation of the approximated surfaces from base points: %.2f and %.2f' %(diff1_range,diff2_range))

        if doPlot:
            plot_regression(pts_1,pts_2,box,poly_reg=[poly_reg_1,poly_reg_2])

        return [box, spec_limits, poly_reg]

    elif mode==SPLINE_APPR_MODE:    # Spline Regression of the boundary surfaces

        log.info('Spline Approximation for the boundary surfaces!')

        spline_reg1 = scipy.interpolate.bisplrep(pts_1[:,0],pts_1[:,1],pts_1[:,2],s=2)
        spline_reg2 = scipy.interpolate.bisplrep(pts_2[:,0],pts_2[:,1],pts_2[:,2],s=2)
        spline_reg = [spline_reg1,spline_reg2]

        if doPlot:
            plot_regression(pts_1,pts_2,box,spline_reg=[spline_reg1,spline_reg2])

        return [box, spec_limits, spline_reg]

    elif mode==INTERP_APPR_MODE:

        log.info('Interpolation method to approach surfaces!')
        
        shape = (X_INTERP,Y_INTERP)
        
        pts_1[:,0] *= float(X_INTERP)/xmax
        pts_2[:,0] *= float(X_INTERP)/xmax
        
        pts_1[:,1] *= float(Y_INTERP)/ymax
        pts_2[:,1] *= float(Y_INTERP)/ymax
        
        order = numpy.asarray(order)

        log.info('Polynomial Regression with order %s!' %(order))

        parameters_1,base_1 = pts_1[:,:2],pts_1[:,2]
        poly_reg_1 = util.polynomial_regression( parameters_1, base_1, order )

        parameters_2,base_2 = pts_2[:,:2],pts_2[:,2]
        poly_reg_2 = util.polynomial_regression( parameters_2, base_2, order )

        bg1 = poly_reg_1.eval(parameters_1)
        bg2 = poly_reg_2.eval(parameters_2)
        
        pts_1[:,2] -= bg1
        pts_2[:,2] -= bg2

        cache_name = trk_path + '.surf'
        cachefile = mrc.MRCFile(cache_name)

        if os.path.exists(cache_name):
            log.info('Initialize interpolation maps from %s' %cache_name)
            init1 = cachefile.getZSliceAt(0)
            init2 = cachefile.getZSliceAt(1)
            map1 = util.mapSurface(shape,pts_1,init=init1)
            map2 = util.mapSurface(shape,pts_2,init=init2)
        else:
            map1 = util.mapSurface(shape,pts_1)
            map2 = util.mapSurface(shape,pts_2)
            cachefile.setHeader(X_INTERP,Y_INTERP,2)

        log.info('Save interpolation maps into %s' %cache_name)
        cachefile.setZSliceAt(0,map1)
        cachefile.setZSliceAt(1,map2)
        cachefile.updateHeader()
                
        # Evaluate the background for the map
        
        x,y = numpy.mgrid[1:X_INTERP+1,1:Y_INTERP+1]
        xy = numpy.column_stack((x.ravel(),y.ravel()))
        
        bg1_map = poly_reg_1.eval(xy)
        bg2_map = poly_reg_2.eval(xy)
        
        bg1_map.resize((X_INTERP,Y_INTERP))
        bg2_map.resize((X_INTERP,Y_INTERP))
        
        # Add the background values back...
        
        pts_1[:,2] += bg1
        pts_2[:,2] += bg2        
        
        map1 += bg1_map
        map2 += bg2_map
        
        # Resize back the map
        
        pts_1[:,0] *= float(xmax)/X_INTERP
        pts_2[:,0] *= float(xmax)/X_INTERP
        
        pts_1[:,1] *= float(ymax)/Y_INTERP
        pts_2[:,1] *= float(ymax)/Y_INTERP
        
        map1_ = util.Resize(map1,(xmax,ymax))
        map2_ = util.Resize(map2,(xmax,ymax))
        
#        map1_ = numpy.zeros((xmax,ymax),dtype='float32')
#        map2_ = numpy.zeros((xmax,ymax),dtype='float32')
#        
#        map1_ = util.array2cv(map1_)
#        map2_ = util.array2cv(map2_)
#        
#        cv.Resize(util.array2cv(map1.astype('float32')),map1_,cv.CV_INTER_CUBIC)
#        cv.Resize(util.array2cv(map2.astype('float32')),map2_,cv.CV_INTER_CUBIC)
#        
#        map1_ = util.cv2array(map1_)
#        map2_ = util.cv2array(map2_)
#        
#        map1_.resize((xmax,ymax))
#        map2_.resize((xmax,ymax))
        
        map = (map1_,map2_)

        if doPlot:
            
            # 3D plots
            
            pts = numpy.row_stack((pts_1,pts_2))

            step = 10
            x,y = numpy.mgrid[1:xmax:step,1:ymax:step]
            
            surf1 = (x,y,map1_[::step,::step].copy())
            surf2 = (x,y,map2_[::step,::step].copy())
            
            s = (float(xmax),float(ymax),float(zmax))
                        
            txbr.utilities.plot3DMarkers(pts[:,0],pts[:,1],pts[:,2],scale=s,surfaces=(surf1,surf2))
            
            # Contour plots
            
            util.plot.plot_contour(map1,title='Contour \#0')
            util.plot.plot_contour(map2,title='Contour \#1')
            
            pylab.show()
            

        return ( box, spec_limits, map)

    return None


def eval_width(limits,frame,regression,mode,width_app=ISO_VOLUME_WIDTH):
    '''This routine evaluate the width of the flattened specimen from the boundary
    approximations of the warped sample (specified by variable mode). Setting the
    width corresponds to setting the new specimen limits. Three type of evaluations
    can be used. If width_app=ISO_VOLUME_WIDTH, volume is kept constant between the
    warped and the flattened specimen. If width_app=MINIMUM_WIDTH, the width of the
    flattened volume in the Z direction corresponds to the width of the initial
    volume in the same direction. Similar case for width_app=MAXIMUM_WIDTH, but with
    a maximum approximation.
    '''

    log.info('Evaluate the final specimen width.')

    zlim_min,zlim_max = limits[2]
    full_wdth = zlim_max-zlim_min

    speclim = lambda width : [int(full_wdth/2.0-width/2.0),int(full_wdth/2.0+width/2.0)]

    xmin,xmax = frame[0]
    ymin,ymax = frame[1]

    # Trim again fin 10%

    delta_x = xmax-xmin
    delta_y = ymax-ymin

    xmin = int(xmin + 0.05*delta_x)
    xmax = int(xmax - 0.05*delta_x)
    ymin = int(ymin + 0.05*delta_y)
    ymax = int(ymax - 0.05*delta_y)

    # Let's do it

    delta = {}

    if mode==POL_APPR_MODE:

        poly_reg_1,poly_reg_2 = regression

        x = numpy.linspace(xmin,xmax,xmax-xmin+1)
        y = numpy.zeros((xmax-xmin+1))

        delta_z_min = numpy.inf
        delta_z_max = -numpy.inf
        delta_z_isoV = 0.0

        for iy in range(ymin,ymax+1):
            y[:] = iy
            z1_reg = poly_reg_1.eval(numpy.column_stack((x,y)))
            z2_reg = poly_reg_2.eval(numpy.column_stack((x,y)))
            delta_z_min = min(delta_z_min,numpy.min(z2_reg-z1_reg))
            delta_z_max = max(delta_z_max,numpy.max(z2_reg-z1_reg))
            delta_z_isoV += numpy.sum(z2_reg-z1_reg)

        delta[ISO_VOLUME_WIDTH] = delta_z_isoV/(xmax-xmin+1)/(ymax-ymin+1)
        delta[MAXIMUM_WIDTH] = delta_z_max
        delta[MINIMUM_WIDTH] = delta_z_min

    if mode==SPLINE_APPR_MODE:

        sp1,sp2 = regression

        x = numpy.linspace(xmin,xmax,xmax-xmin+1)
        y = numpy.zeros((xmax-xmin+1))

        delta_z_min = numpy.inf
        delta_z_max = -numpy.inf
        delta_z_isoV = 0.0

        for iy in range(ymin,ymax+1):
            y[:] = iy
            z1_reg = scipy.interpolate.bisplev(x,y,sp1)
            z2_reg = scipy.interpolate.bisplev(x,y,sp2)
            delta_z_min = min(delta_z_min,numpy.min(z2_reg-z1_reg))
            delta_z_max = max(delta_z_max,numpy.max(z2_reg-z1_reg))
            delta_z_isoV += numpy.sum(z2_reg-z1_reg)

        delta[ISO_VOLUME_WIDTH] = delta_z_isoV/(xmax-xmin+1)/(ymax-ymin+1)
        delta[MAXIMUM_WIDTH] = delta_z_max
        delta[MINIMUM_WIDTH] = delta_z_min

    if mode==INTERP_APPR_MODE:

        map1,map2 = regression

        x = numpy.linspace(xmin,xmax,xmax-xmin+1)
        y = numpy.zeros((xmax-xmin+1))

        delta_z_min = numpy.inf
        delta_z_max = -numpy.inf
        delta_z_isoV = 0.0

        for iy in range(ymin,ymax+1):
            y[:] = iy
            z1_reg = map1[xmin:xmax,iy]
            z2_reg = map2[xmin:xmax,iy]
            delta_z_min = min(delta_z_min,numpy.min(z2_reg-z1_reg))
            delta_z_max = max(delta_z_max,numpy.max(z2_reg-z1_reg))
            delta_z_isoV += numpy.sum(z2_reg-z1_reg)

        delta[ISO_VOLUME_WIDTH] = delta_z_isoV/(xmax-xmin+1)/(ymax-ymin+1)
        delta[MAXIMUM_WIDTH] = delta_z_max
        delta[MINIMUM_WIDTH] = delta_z_min

    log.info('width:[%f,%f,%f]' %(delta[MINIMUM_WIDTH],delta[ISO_VOLUME_WIDTH],delta[MAXIMUM_WIDTH]))

    width = delta[width_app]

    return speclim(width),width


def do_warp(vol_path,regression,spec_limits,mode='polynom',direction='Z',test=False):
    '''Given the file path to the volume to flatten, this routine uses the polynomial or spline
    regression to unwarp the volume'''

    if vol_path==None or not os.path.exists(vol_path):
        raise IOError, 'File %s does not exists!' %vol_path

    # First Calculate the volume between the two boundaries

    Zmin,Zmax = spec_limits

    log.info('Flattening will occur between slices #%i-#%i' %(Zmin,Zmax))

    if mode=='polynom':
        pz1 = regression[0]
        pz2 = regression[1]
    elif mode=='spline':
        sp1 = regression[0]
        sp2 = regression[1]
    elif mode=='interpolation':
        map1 = regression[0]
        map2 = regression[1]

    file = mrc.MRCFile(vol_path)
    file_out = mrc.MRCFile(vol_path + '.uwrpd',template=vol_path)

    if direction=='X':
        n1 = file.ny
        n2 = file.nx
        n3 = file.nz
        slicaxis = 'Z'
    elif direction=='Y':
        n1 = file.nz
        n2 = file.ny
        n3 = file.nx
        slicaxis = 'X'
    else:
        n1 = file.nx
        n2 = file.nz
        n3 = file.ny
        slicaxis = 'Y'

    slice = numpy.zeros((n1,2))
    slice[:,0] = numpy.arange(1.0,float(n1)+1.0)

    z_inc = numpy.arange(n2)

    if test:
        iy0 = int(n3/2.0)
        y_range = range(iy0,iy0+1)
    else:
        y_range = range(n3)

    counter = 0

    index = numpy.zeros((n1,n2),dtype='int')
    delta = numpy.zeros((n1,n2),dtype='float')

    for iy in y_range:

        t0 = time.clock()

        slice[:,1] = iy

        if mode==POL_APPR_MODE:
            Z1 = pz1.eval(slice)
            Z2 = pz2.eval(slice)
        elif mode==SPLINE_APPR_MODE:
            Z1 = scipy.interpolate.bisplev(slice[:,0],[iy],sp1)
            Z2 = scipy.interpolate.bisplev(slice[:,0],[iy],sp2)
        elif mode==INTERP_APPR_MODE:
            Z1 = map1[:,iy]
            Z2 = map2[:,iy]

        Z1 = Z1.ravel()
        Z2 = Z2.ravel()

        indexOfZ1 = numpy.asarray(Z1)
        indexOfZ2 = numpy.asarray(Z2)

        for i2 in range(n2):
            index_ = (indexOfZ2[:]-indexOfZ1[:])*(i2-Zmin)/(Zmax-Zmin) + indexOfZ1[:]
            index[:,i2] = index_    # cast to the floor integer
            delta[:,i2] = index_ - index[:,i2]

        index = numpy.where(index>=0,index,0)
        index = numpy.where(index<n2-1,index,0)

        if direction=='X':
            u = file.getZSliceAt(iy)
            u = u.swapaxes(0,1)
        elif direction=='Y':
            u = file.getXSliceAt(iy)
            u = u.swapaxes(0,1)
        else:
            u = file.getYSliceAt(iy)

        v = numpy.empty_like(u)

        # Make sure interpolation during warping will be the same for every channel (in case of RGB)

        n_channels = u.size/n1/n2
        delta_ = delta.repeat([n_channels]*n1*n2)
        delta_.resize(u.shape)

        for i1 in range(n1):
            zrange = index[i1,:]
            v[i1,:] = (1.0-delta_[i1,zrange])*u[i1,zrange] + delta_[i1,zrange]*u[i1,zrange+1]

        if direction=='X':
            file_out.setZSliceAt(iy,v.swapaxes(0,1))
        elif direction=='Y':
            file_out.setXSliceAt(iy,v.swapaxes(0,1))
        else:
            file_out.setYSliceAt(iy,v)

        t1 = time.clock()

        counter = counter + 1

        delta_t = t1-t0

        log.info('Slice #%i/%i at %s=%i flattened (%fs)' %(counter,len(y_range),slicaxis,iy,delta_t))

        if test:

            title = 'Warped and Flattened Section at %s=%i.' %(slicaxis,iy0)
            plot_sections(u,v,[Z1,Z2],title,mode=file.mode)


def flatten(vol_pth, trk_path, direction='Z', mode='polynom', order=None, width_app="isov", test=False, doPlot=False):
    '''The routine implements a volume flattening. The MRC volume file is specified in the file system
    by vol_pth). Boundaries
    need to be egmented with points in a model file, one object for each of the two surfaces. Points
    can be in different contours. Those surfaces are
    then approximated by polynomial functions, or with splines or interpolated between the input data
    points according to the 'mode' variable ('polynom','spline','interpolation').
    mode: (polynom,spline,interpolation)
    '''

    log.info('Execute sample flattening in the %s direction' %direction)
    
    if order==None and mode=='polynom': order = [10,10]
    if order==None and mode=='interpolation': order = [2,2]

    keywords = {'direction':direction,'mode':mode,'order':order,'doPlot':doPlot}
    (limits,frame,regression) = eval_spl_boundaries(trk_path,**keywords)    # z deviation

    (speclim,width) = eval_width(limits,frame,regression,mode,width_app=width_app)

    keywords = {'direction':direction,'mode':mode,'test':test}
    do_warp(vol_pth,regression,speclim,**keywords)


def plot_regression(pts1,pts2,limits,poly_reg=None,spline_reg=None):
    '''Plot regression of the specimen boundary surfaces'''

    xmin,xmax = limits[0]
    ymin,ymax = limits[1]

    from enthought.mayavi import mlab

    fig = mlab.figure()

    mlab.points3d(pts1[:,0],pts1[:,1],pts1[:,2],scale_factor=20.0)
    mlab.points3d(pts2[:,0],pts2[:,1],pts2[:,2],scale_factor=20.0)

    if poly_reg!=None:

        poly_reg_1,poly_reg_2 = poly_reg

        z1_reg = poly_reg_1.eval(numpy.column_stack((pts1[:,0],pts1[:,1])))
        z2_reg = poly_reg_2.eval(numpy.column_stack((pts2[:,0],pts2[:,1])))

        mlab.points3d(pts1[:,0],pts1[:,1],z1_reg,scale_factor=20.0,color=(0.0,1.0,1.0))
        mlab.points3d(pts2[:,0],pts2[:,1],z2_reg,scale_factor=20.0,color=(0.0,1.0,1.0))

    if spline_reg!=None:

        spline_reg1,spline_reg2 = spline_reg

        xnew,ynew = scipy.mgrid[100:xmax-100:50j,100:ymax-100:50j]

        z1_spline = scipy.interpolate.bisplev(xnew[:,0],ynew[0,:],spline_reg1)
        z2_spline = scipy.interpolate.bisplev(xnew[:,0],ynew[0,:],spline_reg2)

        xnew = xnew.ravel()
        ynew = ynew.ravel()

        z1_spline = z1_spline.ravel()
        z2_spline = z2_spline.ravel()

        mlab.points3d(xnew,ynew,z1_spline,scale_factor=5.0,color=(0.0,1.0,0.0))
        mlab.points3d(xnew,ynew,z2_spline,scale_factor=5.0,color=(0.0,1.0,0.0))

    mlab.show()


def plot_sections(u,v,limits,title,mode=None):

        u = u.T
        v = v.T
        (Z1,Z2) = limits

        (n1,n2) = u.shape

        kws = {}
        if mode!=16:
            kws['cmap'] = matplotlib.cm.gray

        pylab.figure()

        pylab.subplot(211)
        pylab.imshow(u,**kws)
        pylab.plot(Z1)
        pylab.plot(Z2)
        pylab.xlim(0,n1)
        pylab.ylim(0,n2)

        pylab.subplot(212)
        pylab.imshow(v,**kws)
        pylab.xlim(0,n1)
        pylab.ylim(0,n2)

        pylab.suptitle(title)

        pylab.show()
