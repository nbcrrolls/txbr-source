import os.path
import os
import numpy
import numpy.linalg
import numpy.fft
import scipy.signal
import pylab
import mrc
import re
import util
import scipy.ndimage.measurements
import scipy.ndimage.filters

numpy.set_printoptions(precision=3)


def changeToLogAxis(u):

    nx,ny = u.shape

    x = numpy.logspace(1,numpy.log10(nx),nx) - 1.0
    y = numpy.logspace(1,numpy.log10(ny),ny) - 1.0

    x0 = x.astype('int')
    y0 = y.astype('int')

    px = x - x0
    py = y - y0

    v = numpy.empty_like(u)

    for index,x_ in enumerate(x0[:-1]):
        v[index,:] = px[index]*u[x_+1,:] + (1.0-px[index])*u[x_,:]
    v[-1,:] = u[-1,:]

    u = v.copy()

    for index,y_ in enumerate(y0[:-1]):
        v[:,index] = py[index]*u[:,y_+1] + (1.0-py[index])*u[:,y_]
    v[:,-1] = u[:,-1]

    return v


def maximum( u , mode="gaussian-fit", verbose=False, doPlot=False):
    '''Characterizes the ("correlation") peak in an image array
    Input parameters:
        u: 2D numpy array containing the data
        mode: in ["gaussian-fit","parabolic"]
    '''

    u = numpy.asarray(u).copy()
    u = numpy.fft.fftshift(u)  # To bring the maximum roughly in the center of the image
    u = u.real

    nx,ny = u.shape

    max_loc = numpy.argmax(u)   # for the flattened array

    max_row = max_loc/ny
    max_col = max_loc - max_row*ny
    peak = numpy.max(u)

    tx = max_row - nx/2
    ty = max_col - ny/2

    if verbose:
        print "Maximum: %s at (%i,%i) on a %ix%i domain" %(str(numpy.max(u)),max_row,max_col,nx,ny)

    L = 15

    if nx>2*L and ny>2*L:
        centerx,centery = max_row,max_col
        sx,sy = L,L
    else:
        centerx,centery = nx/2,ny/2
        sx,sy = nx/2,ny/2

    #print "Max=(%i,%i) (nx,ny)=(%.1f,%.1f)  (sx,sy)=(%.1f,%.1f)" %(max_row,max_col,nx,ny,sx,sy)

    x1, x2 = max(0, centerx-sx), min(nx, centerx + sx)
    y1, y2 = max(0, centery-sy), min(ny, centery + sy)

    data = u[x1:x2,y1:y2]
    data = data - numpy.min(data)

    y,x = numpy.meshgrid(range(y2-y1),range(x2-x1))
    r2 = (x-(x2-x1)/2.0)**2 + (y-(y2-y1)/2.0)**2

    import txbr.onthefly.otfutil

    alpha = 1.0/2.0/numpy.max((x2-x1,y2-y1))
    std = numpy.sum(data*numpy.exp(-alpha*r2))/numpy.sum(data)

    std = scipy.ndimage.measurements.standard_deviation(data)

    if mode=="parabolic":

        # parabolic interpolation

        dx = (u[max_row+1,max_col] - u[max_row-1,max_col])/2.0
        dy = (u[max_row,max_col+1] - u[max_row,max_col-1])/2.0

        dx2 = u[max_row+2,max_col] - 2*u[max_row,max_col] + u[max_row-1,max_col]
        dy2 = u[max_row,max_col+2] - 2*u[max_row,max_col] + u[max_row,max_col-1]
        dxdy = u[max_row+1,max_col+1] - 2*u[max_row,max_col] + u[max_row-1,max_col-1]

        A = numpy.array([[dx2,dxdy],[dxdy,dy2]])
        B = -numpy.array([dx,dy])

        p = numpy.dot(numpy.linalg.inv(A),B)

        px = p[0]
        py = p[1]

        # Translation

        tx = max_row + px - nx/2
        ty = max_col + py - ny/2

        if verbose:
            print "Parabolic approximation for the peak"
            print "Maximum: (%i,%i)  (%.1f,%.1f)" %(max_row,max_col,tx,ty)
            print "dx2=%f   dy2=%f   dxdy=%f" %(dx2,dy2,dxdy)

        return (tx,ty)

    elif mode=="gaussian-fit":
        
#        L = 25
#
#        if nx>2*L and ny>2*L:
#            centerx,centery = max_row,max_col
#            sx,sy = L,L
#        else:
#            centerx,centery = nx/2,ny/2
#            sx,sy = nx/2,ny/2
#
#        #print "Max=(%i,%i) (nx,ny)=(%.1f,%.1f)  (sx,sy)=(%.1f,%.1f)" %(max_row,max_col,nx,ny,sx,sy)
#
#        x1, x2 = max(0, centerx-sx), min(nx, centerx + sx)
#        y1, y2 = max(0, centery-sy), min(ny, centery + sy)
#
#        data = u[x1:x2,y1:y2]
        
        data = data.T

        params = [numpy.max(data),0.0,200.0,200.0,2.0,2.0,90.0]
        usemoment = [True,True,True,True,True,True,False]

        try:
            params = util.gaussfit(data,params=params,usemoment=usemoment)
        except:
            print "Max=(%i,%i) (nx,ny)=(%.1f,%.1f)  (sx,sy)=(%.1f,%.1f)" %(max_row,max_col,nx,ny,sx,sy)
            return 0.0,0.0,1.0,1.0,90.0

        b,a,x,y,width_x,width_y,angle = params

        if width_x>width_y:
            width_x,width_y = width_y,width_x
            angle = angle + 90
            params = b,a,x,y,width_x,width_y,angle

        if numpy.abs(x)>2*sx or numpy.abs(y)>2*sy:
            print "Out of range problem with Maximum!"
            x = sx
            y = sy
            a = 0.0 # Notify the peak was not found

        angle = (angle+90)%180-90

        if verbose:
            print "Gaussan-Fit: (a,b)=(%.1e,%.1e)  center=(%.2f,%.2f)  width=(%.2f,%.2f)" %(a,b,x,y,width_x,width_y)
            print "             Angle in degrees: %.2f  [cos=%.2f,sin=%.2f]" %(angle,numpy.cos(numpy.radians(angle)),numpy.sin(numpy.radians(angle)))

        tx = centerx - sx + x - nx/2
        ty = centery - sy + y - ny/2

        if doPlot:

            pylab.matshow(data, cmap=pylab.cm.gist_earth_r)
            fit = util.twodgaussian(params)
            pylab.contour(fit(*pylab.indices(data.shape)), cmap=pylab.cm.copper)
            ax = pylab.gca()
            pylab.show()

        #return ( tx, ty, width_x, width_y, angle )

        return ( tx, ty, width_x, width_y, angle, a, b )

    return ( tx, ty, std, peak )


def getPhaseCorrelation( u1, u2, apply_med_filt=0 ):
    '''Create the phase correlation map between two arrays u1 and u2'''

    apply_hamming = True    # To remove artefacts in k=0
    apply_gauss = True
    apply_exp = True
    apply_median_filter = False
    apply_median_filter = apply_med_filt!=0

    if len(u1.shape)==1:

        nx = u1.size

        if apply_hamming:
            u1 = u1*numpy.hamming(nx)
            u2 = u2*numpy.hamming(nx)

        u1 = numpy.fft.fft(u1)
        u2 = numpy.fft.fft(u2)

        u =  u1*numpy.conj(u2)/numpy.abs( u1*u2 )

        u = numpy.fft.ifft(u)

        if apply_median_filter: u = scipy.signal.medfilt(numpy.abs(u),kernel_size=apply_med_filt)
        if apply_exp: u = numpy.exp(numpy.abs(u))-1.0

    if len(u1.shape)==2:

        nx,ny = u1.shape

        if apply_hamming:
            d = numpy.outer(numpy.hamming(nx),numpy.hamming(ny))
            u1 = u1*d
            u2 = u2*d

        u1 = numpy.fft.fft2(u1)
        u2 = numpy.fft.fft2(u2)

        u =  u1*numpy.conj(u2)/numpy.abs( u1*u2 )

        if apply_gauss:
            std = 100.0
            gauss = numpy.outer(scipy.signal.gaussian(nx,std),scipy.signal.gaussian(ny,std))
            u = u*numpy.fft.ifftshift(gauss)

        u = numpy.fft.ifft2(u)

        if apply_exp: u = numpy.exp(numpy.abs(u.real))-1.0
        if apply_median_filter: u = scipy.signal.medfilt2d(numpy.abs(u),kernel_size=apply_med_filt)

    return u


def correlate( u1, u2, crop_rate=0.0, mode='gaussian-fit', fullInfo=False, verbose=False ):

    if crop_rate!=0.0:

        x_margin = min(int(crop_rate*u1.shape[0]),int(crop_rate*u2.shape[0]))
        y_margin = min(int(crop_rate*u1.shape[1]),int(crop_rate*u2.shape[1]))

        x_margin = max(1,x_margin)
        y_margin = max(1,y_margin)

        u1_ = u1[x_margin:-x_margin,y_margin:-y_margin]
        u2_ = u2[x_margin:-x_margin,y_margin:-y_margin]

    else:

        u1_ = u1
        u2_ = u2


    phase = getPhaseCorrelation( u1_, u2_ )
    
    peak_loc = maximum( phase, mode=mode, verbose=verbose, doPlot=False )

    if fullInfo:
        return peak_loc,phase
    else:
        return peak_loc
 

def prealignFile( file, tilts=None, crop_rate=0.0, createAlignedStack=True, createCorrelationStack=False ):
    '''Prealign an MRC file using phase correlation.

    Arguments:
    "file": The MRC file to pre-align
    "tilts": A list of the index to consider in the original file. If set to None
             every slices will be taken into account
    index -- Index of the under-focused image within the MRC file
    '''

    if not os.path.lexists(file):
        print "File %s does not exist!" %file
        return;

    base_match = re.match("(\S*)\\.(\w*)$", file)
    if base_match:
        basename = os.path.basename(base_match.group(1))
    else:
        basename = None

    f = mrc.MRCFile(file)

    if tilts==None:
        tilts = range(f.nz)
    else:
        tilts = [ tilt for tilt in tilts if tilt<f.nz ]

    x_margin = max(1,int(crop_rate*f.nx/2.0))
    y_margin = max(1,int(crop_rate*f.ny/2.0))

    corr_info = []
    corr_shifts = []

    for index,tilt1 in enumerate(tilts):
        if index==len(tilts)-1: continue
        tilt2 = tilts[index+1]
        if (tilt1>=f.nz-1): break
        u1 = f.getZSliceAt(tilt1)
        u2 = f.getZSliceAt(tilt2)
        u1 = u1[x_margin:-x_margin,y_margin:-y_margin]
        u2 = u2[x_margin:-x_margin,y_margin:-y_margin]
#        u1 = changeToLogAxis(u1)
#        u2 = changeToLogAxis(u2)
        phase_corr = getPhaseCorrelation( u1, u2 )
        corr_info.append(maximum(phase_corr,doPlot=False))
        corr_shifts.append(numpy.fft.fftshift(phase_corr))
        print "Tilt sequence %3i/%-3i: t=(%+6.3f,%+6.3f)  width=(%+7.3f,%+7.3f) ratio:%.3f  angle=%-3.3f" %(tilt1,tilt2,corr_info[-1][0],corr_info[-1][1],corr_info[-1][2],corr_info[-1][3],corr_info[-1][3]/corr_info[-1][2],corr_info[-1][4])

    corr_info = numpy.asarray(corr_info)

    # Save the correlation info

    numpy.savetxt( "%s_corr_info.txt" %basename, corr_info )

    # PLot the angles

    pylab.figure()
    pylab.plot(corr_info[:,4])
    pylab.show()

    # Create the corresponfing transformations

    translations = numpy.vstack(([[0.0,0.0]],corr_info[:,:2]))
    translations = numpy.cumsum(translations,axis=0)

    index_ref = translations.shape[0]/2
    tr_ref = translations[index_ref]

    translations = translations - tr_ref

    transformations = {}
    for index, tilt1 in enumerate(tilts):
        transformations[tilt1] = util.Translation2D(*translations[index])

    # End of creating transformations

    if createCorrelationStack: # Eventually create a correlation stack...
        mrc.createMRCFile( "%s_corr.mrc" %basename, numpy.dstack(corr_shifts) )
        os.system("imod %s_corr.mrc" %basename)

    if createAlignedStack: # Eventually create an aligned stack...
        f.warp( transformations, output="%s.preali" %basename)
        os.system("imod %s.preali" %basename)


def test():

    nx = 500
    ny = 500

    alpha1 =  numpy.array([ 0.1, 0.1 ])
    
    x1 = numpy.arange(-nx/2.0,nx/2.0,1.0)
    y1 = numpy.arange(-ny/2.0,ny/2.0,1.0)

    u1 = numpy.outer(numpy.exp(-alpha1[0]*x1**2),numpy.exp(-alpha1[1]*y1**2))

    shift2 = numpy.array([ 0.0, 0.0 ])
    alpha2 = numpy.array([ 0.02, 0.02 ])
    alpha2 = 0.5*alpha1

    x2 = numpy.arange(-nx/2.0+shift2[0],nx/2.0+shift2[0],1.0)
    y2 = numpy.arange(-ny/2.0+shift2[1],ny/2.0+shift2[1],1.0)
    
    u2 = numpy.outer(numpy.exp(-alpha2[0]*x2**2),numpy.exp(-alpha2[0]*y2**2))

    u1 = u1/numpy.std(u1)
    u2 = u2/numpy.std(u2)

    u1 = changeToLogAxis(u1)
    u2 = changeToLogAxis(u2)
#
    phase_corr = getPhaseCorrelation( u1, u2 )
 #   print maximum(phase_corr)

    u3 = numpy.fft.fftshift(phase_corr)
    u3 = u3/numpy.std(u3)
    
    mrc.createMRCFile( 'output.mrc', numpy.dstack((u1,u2,u3)).real )

    os.system("imod output.mrc")


if __name__ == '__main__':

    file = "/home/sphan/data/MMV-7/MMV-7.st"
    #file = "/home/sphan/data/MMV-7/fhv6a.st"
    #file = "/home/sphan/data/MMV-7/FHV102709-29.st"
    #file = "/home/sphan/data/MMV-7/MMV-9.st"
    file= "/ncmirdata3/mterada/for_sebs/2Dbin8.st"
    file= "/ncmirdata3/sphan/otf/14otf/test.st"
    #file= "/ncmirdata3/sphan/otf/14otftiezt/tietz.st"
    file= "/ncmirdata3/sphan/otf/15de12/FHV-16_de12.st"
    file= "/ncmirdata3/sphan/otf/15de12/FHV-16_de12.st"
    file= "/ncmirdata3/sphan/FHV-18/bin8/FHV-18a.st"


    prealignFile( file, createAlignedStack=True, createCorrelationStack=True )

    #test()

