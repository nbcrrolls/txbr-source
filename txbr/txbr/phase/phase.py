import os.path
import os
import sys
import numpy
import numpy.fft
import cv
import mrc
import modl
import util

global ROW_MARGIN, COL_MARGIN

ROW_MARGIN = 500
COL_MARGIN = 500

def createMRCFile( filename, u ):
    '''Create an MRC file from a 3D numpy array'''

    nx,ny,nz = u.shape

    print "MRC File: (nx,ny,nz)=(%i,%i,%i)" %(nx,ny,nz)

    output = mrc.MRCFile(filename)
    output.setHeader(nx,ny,nz)
    for iz in range(nz):
        output.setZSliceAt(iz,u[:,:,iz])
    output.updateHeader()

    os.system("imod %s" %filename)


def align( uref, u, ptsref=None, pts=None ):
    '''Align u on uref.

    Arguments:
    uref -- The image reference (as a 2D numpy array)
    u -- The image (as a 2D numpy array) to be aligned

    '''
    
    doCrossCorrelation = ptsref==None or pts==None or len(ptsref)!=len(pts) or len(ptsref)<3

    if doCrossCorrelation:
        return alignFromCrossCorrelation( uref, u )
    else:
        return alignFromPoints( uref, u, ptsref, pts)
        

def alignFromPoints( uref, u, ptsref, pts):
    '''Align u on uref by translating it to match the largest cross-correlation peak.
    Return a copy of the translated u.

    Arguments:
    uref -- The image reference (as a 2D numpy array)
    u -- The image (as a 2D numpy array) to be aligned
    ptsref -- Points of the reference image (as a 2D numpy array)
    pts -- Corresponding points in the image to align

    '''
    
    print "Starting Image alignment by using corresponding points."

    mat1 = cv.fromarray(uref.astype('float32'))
    mat2 = cv.fromarray(u.astype('float32'))

    nx,ny = cv.GetDims(mat1)

    # Create the warping map

    useOrthogonalMap = False

    if useOrthogonalMap:
        R, theta, tx, ty = util.orthogonal_regression_2d( pts, ptsref )
        coeff = R.coeff
    else:
        P1, P2 = util.polynomial_regression( pts, ptsref, [1,1] )
        coeff = []
        coeff.extend(P1.coeff)
        coeff.extend(P2.coeff)

    warp = cv.CreateMat(2,3,cv.CV_32FC1)
    # Axes are shifted between openCV and numpy

    warp[0,2] = coeff[3]
    warp[1,2] = coeff[0]

    warp[0,0] = coeff[5]
    warp[0,1] = coeff[4]
    warp[1,0] = coeff[2]
    warp[1,1] = coeff[1]

    print "warp x: %10.4f    %10.4f %10.4f" %(warp[0,2],warp[0,0],warp[0,1])
    print "warp y: %10.4f    %10.4f %10.4f" %(warp[1,2],warp[1,0],warp[1,1])

    # Apply the warping map

    dest = cv.CreateMat(nx,ny,cv.GetElemType(mat2))
    cv.WarpAffine(mat2, dest, warp)

    print "End Image alignment."

    return numpy.asarray(dest,dtype='float')

    
def alignFromCrossCorrelation( uref, u ):
    '''Align u on uref by translating it to match the largest cross-correlation peak.
    Return a copy of the translated u.

    Arguments:
    uref -- The image reference (as a 2D numpy array)
    u -- The image (as a 2D numpy array) to be aligned

    '''

    print "Starting Image alignment by cross-correlation."

    mat1 = cv.fromarray(uref.astype('float32'))
    mat2 = cv.fromarray(u.astype('float32'))

    print "Margin: (mx,my)=(%i,%i)" %(ROW_MARGIN,COL_MARGIN)

    nx,ny = cv.GetDims(mat1)

    mat2_ = cv.GetSubRect(mat2,(ROW_MARGIN, COL_MARGIN, mat2.rows-2*ROW_MARGIN, mat2.cols-2*COL_MARGIN))
    template = cv.CreateMat(2*ROW_MARGIN+1, 2*ROW_MARGIN+1, cv.CV_32FC1)

    cv.MatchTemplate( mat1, mat2_, template, cv.CV_TM_CCORR_NORMED )

    minVal, maxVal, minLoc, maxLoc = cv.MinMaxLoc(template)

    t = numpy.array(maxLoc) - numpy.array((ROW_MARGIN, COL_MARGIN))

    if t[0]>=ROW_MARGIN or t[1]>=COL_MARGIN:
        print "Margins are too small"
        sys.exit(0)

    # Create the warping map

    warp = cv.CreateMat(2,3,cv.CV_32FC1)

    warp[0,2] = t[0]
    warp[1,2] = t[1]

    warp[0,0] = 1.0
    warp[0,1] = 0.0
    warp[1,0] = 0.0
    warp[1,1] = 1.0

    print "warp x: %.2f    %.2f %.2f" %(warp[0,2],warp[0,0],warp[0,1])
    print "warp y: %.2f    %.2f %.2f" %(warp[1,2],warp[1,0],warp[1,1])

    # Apply the warping map

    dest = cv.CreateMat(nx,ny,cv.GetElemType(mat2))
    cv.WarpAffine(mat2, dest, warp)

    print "End Image alignment."

    return numpy.asarray(dest,dtype='float')


def getPhaseFromFile( file, index1=0, index2=2, indexOfFocus=1, model_file=None, output="output.mrc", eps = 1.0e-4, test_align=False, test_diff=False ):
    '''Calculate the phase from three snapshots (underfocused, overfocused and
    in focus image) contained in an MRC file. The phase image is stored in the MRC file
    output

    Arguments:
    file - The MRC file that contains the three snapshots
    index1 -- Index of the under-focused image within the MRC file
    index2 -- Index of the over-focused-image within the MRC file
    indexOfFocus -- Index of the focused image within the MRC file
    output - The MRC file containing the phase result

    '''

    # Load the Slices

    f = mrc.MRCFile(file)

    u1 = f.getZSliceAt(index1)
    u2 = f.getZSliceAt(index2)
    ufocus = f.getZSliceAt(indexOfFocus)
    
    # Load the model points
    
    pts1, pts2, ptsfocus = [], [], []

    if model_file!=None and os.path.exists(model_file):

        model = modl.Model(model_file)
        model.loadFromFile()

        for contour in model.objects[0].contours:
            pts1.append(contour.points[0,:2])
            pts2.append(contour.points[1,:2])
            ptsfocus.append(contour.points[2,:2])

        pts1 = numpy.vstack(pts1)
        pts2 = numpy.vstack(pts2)
        ptsfocus = numpy.vstack(ptsfocus)

    getPhase(u1,u2,ufocus,pts1=pts1,pts2=pts2,ptsfocus=ptsfocus,output=output,eps=eps,test_align=test_align,test_diff=test_diff)



def getPhase( u1, u2, ufocus, pts1=None, pts2=None, ptsfocus=None, output="output.mrc", eps=1.0e-4, test_align=False, test_diff=False ):
    '''Calculate the phase from three snapshots (an underfocused, an overfocused and
    an in focus image).

    Arguments:
    u1 -- The under-focused image (as a 2D numpy array)
    u2 -- The over-focused-image (as a 2D numpy array)
    ufocus -- The focused images (as a 2D numpy array). Used as a reference for alignment
    output - The MRC file containing the phase result
    eps -- Regularization parameter
    
    '''

    nx,ny = u1.shape

    print "Dimensions: (nx,ny)=(%i,%i)" %(nx,ny)

    fx = numpy.fft.fftfreq(nx)
    fy = numpy.fft.fftfreq(ny)

    [fx,fy] = numpy.meshgrid(fx,fy)
    f2 = fx**2 + fy**2

    # Take the difference between the two micrographs

    u1 = align(ufocus,u1,ptsref=ptsfocus,pts=pts1)
    u2 = align(ufocus,u2,ptsref=ptsfocus,pts=pts2)

    onlyPlotAlignedSection = test_align

    if onlyPlotAlignedSection:
        createMRCFile( "phase-align.mrc", numpy.dstack((u1,u2,ufocus)) )
        sys.exit(0)

    u = u2 - u1

    #createMRCFile( "test.mrc", numpy.dstack((u,)) )

    # Zeros the difference

    zeroDiff = True

    if zeroDiff:
        average1 = numpy.average(u)
        u = u - numpy.average(u)
        average2 = numpy.average(u)
        print "Difference average: %.2f (before) -> %.2f (after)" %(average1,average2)

    onlyPlotDiffSection = test_diff

    if onlyPlotDiffSection:
        createMRCFile( "phase-align.mrc", numpy.dstack((u,)) )
        sys.exit(0)

    shortcut = False

    if not shortcut:

        u = numpy.fft.fft2(u)
        u = u/(eps+f2)
        
        ux = fx*u
        uy = fy*u

        ux = numpy.fft.ifft2(ux)
        uy = numpy.fft.ifft2(uy)

        #img = (u1+u2)/2.0   # Approximation of the image

        ux = ux/ufocus
        uy = uy/ufocus

        ux = numpy.fft.fft2(ux)
        uy = numpy.fft.fft2(uy)

        u = (fx*ux+fy*uy)/(eps+f2)
        u = numpy.fft.ifft2(u)

    else:

        u = numpy.fft.fft2(u)
        u = u/(eps+f2)
        u = numpy.fft.ifft2(u)

    print "Average(u): %.2f" %(numpy.average(u))

    # Change the mean and standard deviation of the focused image so both
    # the phase and attenuation image can be displayed at the same time
    
    m_ = numpy.mean(u.real)
    std_ = numpy.std(u.real)

    ufocus = m_ + (ufocus-numpy.mean(ufocus))*std_/numpy.std(ufocus)

    s = 200
    createMRCFile( output, numpy.dstack((u.real[s:-s,s:-s],ufocus[s:-s,s:-s])) )



if __name__ == '__main__':

    file = "/ncmirdata5/sphan/ddd/carbon-grid/defocus_series_3200_30k_092710.mrc"
    indexOfFocus = 11
    index1 = indexOfFocus - 1
    index2 = indexOfFocus + 1

    getPhaseFromFile( file, index1=index1, index2=index2, indexOfFocus=indexOfFocus )

#    import scipy
#    file = "/ncmirdata5/sphan/onthefly2/src_jpegs/HPF3Viewprep070610b_0.0.tif"
#    img = scipy.misc.imread(file)
#    scipy.misc.imshow(img)