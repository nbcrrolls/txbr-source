import os.path
import mrc
import math
import numpy
import scipy
import pylab

from txbr.onthefly.phasecorr import maximum
from txbr.onthefly.phasecorr import getPhaseCorrelation


def _correctForDrift( tranformations, driftOrderCorrection ):
    '''Correct a potential drift in the transformations to align a 3view stack.

    :tranformations: A numpy array of shape (n,3). The first column, contains the
        indices of *good* images, while the second and third column are the actual
        translation vector values to correct for.
    :driftOrderCorrection: The order of the correction to use in the translation variations
    '''

    if len(tranformations)==0:
        return tranformations

    if len(tranformations)<=driftOrderCorrection:
        driftOrderCorrection = len(tranformations)-1

    t = scipy.linspace( 0.0, 1.0, num=len(tranformations) )

    px = scipy.polyfit( t, tranformations[:,1], driftOrderCorrection )
    py = scipy.polyfit( t, tranformations[:,2], driftOrderCorrection )

    tranformations[:,1] -= scipy.polyval( px, t)
    tranformations[:,2] -= scipy.polyval( py, t)

    return tranformations



def alignStack( filename, directory=".", Lx=1000, Ly=1000, driftOrderCorrection=2,
                corr_ratio_threshold = 15.0, doPlotRegistration=False,
                createAlignedStack=True, createBadImageStack=True, verbose=True ):
    '''This function aligns a 3view stack.

    A stack of (3view) images is aligned using pair-wise cross-correlation procedures.
    Bad images are removed using a threshold criteria in the corresponding cross-correlation
    peak between sequential images.
    For performance sake cross-correlation can be performed on a smaller centered frame
    (of half-size *Lx* and *Ly*).
    To remove *drifting* effects, the set of translations can be regressed using a polynomial
    scheme.
    The aligned stack is saved by default into a MRC file having extension ".restack"
    (default usage). Bad images are stored as well in a MRC file having the extension
    ".bad".

    :param filename: Name of the MRC stack to align.
    :param directory: Directory where the MRC stack is located.
    :param Lx, Ly: Half-size of the window frame used for cross-correlation in x and y. If None,
        the whole window is used.
    :param driftOrderCorrection: Order of the polynomial regression used for the drift correction.
        If set to None, no drift correction is performed.
    :param createAlignedStack: Create a final aligned stack, with the bad images excluded.
    :param createBadImageStack: Create a MRC file containing only bad images.

    '''


    filename = os.path.join( directory, filename )

    if not os.path.exists(filename):
        print "The file %s does not exists!" %(filename)
        return
    
    f = mrc.MRCFile(filename)

    if Lx==None: Lx = int(0.25*f.nx)
    if Ly==None: Ly = int(0.25*f.ny)

    slices = range(f.nz)

    x1 = max(0,f.nx/2-Lx)
    x2 = min(f.nx,f.nx/2+Lx)
    y1 = max(0,f.ny/2-Ly)
    y2 = min(f.ny,f.ny/2+Ly)

    badviews = []
    registration = []
    
    for index1,slice in enumerate(slices):
        peak_value = 0.0
        if index1 in badviews: continue
        if index1>=len(slices)-1: break
        index2 = index1
        u1 = f.getZSliceAt(index1)
        survey = True
        while survey and index2<len(slices)-1:
            u2 = f.getZSliceAt(index2+1)
            u1_ = u1[x1:x2,y1:y2]
            u2_ = u2[x1:x2,y1:y2]
            phase = getPhaseCorrelation( u1_, u2_ )
            tx, ty, width_x, width_y, angle, peak_value, background =  maximum( phase, verbose=False, doPlot=False )
            survey = math.fabs(peak_value/background)<corr_ratio_threshold
            if survey:
                badviews.append(index2+1)
                print "Bad View at %i --- %.2f < %.2f " %(index2+1,math.fabs(peak_value/background),corr_ratio_threshold)
            index2 = index2+1
        if not survey: # else we have a problem
            registration.append((index1,index2,tx,ty))
        if verbose:
            print "Peak (%i-%i): %.3e at (%.1f,%.1f) with width (%.1f,%.1f)  [criteria  %.2f >= %.2f]" %( index1, index2, peak_value, tx, ty, width_x, width_y,math.fabs(peak_value/background),corr_ratio_threshold )
        numpy.savetxt("registration.txt",registration)

    print "Number of good views %i/%i" %( len(registration)+1, f.nz )
    print "Bad Views %s" %str(badviews)

    registration = numpy.loadtxt("registration.txt")
    print registration

    # Clean up the registration data to form a set of transformations coordinates

    registration = numpy.asarray(registration)

    if len(registration.shape)==1:
        registration.resize((1,registration.size))

    if registration.shape[1]>1:
        registration[:,2] = numpy.cumsum(registration[:,2])
        registration[:,3] = numpy.cumsum(registration[:,3])

    transformations = numpy.row_stack(([0.0,0.0,0.0],registration[:,1:]))

    # Remove the drift

    driftOrderCorrection = 20
    transformations = _correctForDrift( transformations, driftOrderCorrection )

    doPlotRegistration = True

    if doPlotRegistration:

        print transformations.shape

        pylab.figure()
        pylab.plot(transformations[:,1])
        pylab.plot(transformations[:,2])
        pylab.show()

    # Store data

    if createBadImageStack:

        badslices = numpy.zeros((len(badviews),3))
        badslices[:,0] = numpy.asarray(badviews)

        mrc.translateMRCFile( filename, "%s.bad" %filename, badslices )

    if createAlignedStack:

        print "Creating an aligned stack..."
        
        mrc.translateMRCFile( filename, "%s.restack" %filename, transformations )

    return


if __name__ == '__main__':

    directory = "/ncmirdata3/sphan/"
    filename = "merlin.st"

    options1 = { "directory":"/ncmirdata3/sphan/",  "filename":"merlin.st" }
    options2 = { "directory":"/ncmirdata3/sphan/",  "filename":"merlin_.st" }
    options3 = { "directory":"/ncmirdata3/sphan/",  "filename":"output.mrc" }
    options4 = { "directory":"/ncmirdata3/sphan/",  "filename":"output_.mrc" }

    options = options2

    alignStack( options["filename"], directory=options["directory"] )