import os.path
import txbr
import txbr.utilities
import numpy
import mrc
import cv

def project( p, indexOfSrcSeries, indexOfDestSeries, srcVolume=None, output=None ):

    series = p.series[indexOfDestSeries]
    destFrame = series.getDestinationFrame()

    xmin, ymin, zmin = numpy.min(destFrame,axis=1)
    xmax, ymax, zmax = numpy.max(destFrame,axis=1)
    nx, ny, nz = int(xmax-xmin+1), int(ymax-ymin+1), int(zmax-zmin+1)

    srcSeries = p.series[indexOfSrcSeries]
    fileName =  srcSeries.getReconstructionFileName()

    if not os.path.exists(fileName):
        print "%s does not exist!" %(fileName)
        return

    if srcVolume==None:
        vol = txbr.utilities.MRCStack(fileName)
    else:
        vol = txbr.utilities.MRCStack(srcVolume)

    print "X: %i  %i" %(nx,vol.getImageSize()[0])
    print "Y: %i  %i" %(ny,vol.getImageSize()[1])
    print "Z: %i  %i" %(nz,vol.getNumberOfImages())

    tilts = range(series.numberOfExposures())
    ntilt = len(tilts)


    # ---------------hack ---------

    nx = vol.getImageSize()[0]
    ny = vol.getImageSize()[1]

    # ---------------hack ---------


    if output==None:
        output = "output.mrc"
    output = mrc.MRCFile(output)
    output.setHeader(nx,ny,ntilt)

    for index,itilt in enumerate(tilts):
        u = numpy.zeros((nx,ny))
        map3D = series.projection.getAffineTransformation(itilt)
        print "Tilt #%i: %s" %(itilt,str(map3D))
        for iz in range(nz):
      #      u = u + remap_image(numpy.exp(-vol.getImageAt(iz)),map3D,Z=zmin+iz)
            u = u + remap_image(vol.getImageAt(iz),map3D,Z=zmin+iz)
        output.setZSliceAt(index,u)

    output.updateHeader()

    print "End of the projection routine"


def remap_image( src, map3D, Z=0.0 ):
    '''Helper function to stretch an image (input data: MAT array) that was taken at a given tilt "angle"
    for it to correspond to the slice at Z=0 in the final volume'''
    
    src = src.astype('float32')
    src = cv.fromarray(src)

    nx,ny = cv.GetDims(src)

    dest = cv.CreateMat(nx,ny,cv.GetElemType(src))

    warp = cv.CreateMat(2,3,cv.CV_32FC1)

    # Apply the warp onto the image

    warp[0,0] = map3D.M[1,1]
    warp[0,1] = map3D.M[1,0]
    warp[1,0] = map3D.M[0,1]
    warp[1,1] = map3D.M[0,0]

    warp[0,2] = map3D.T[1] + map3D.M[1,2]*Z
    warp[1,2] = map3D.T[0] + map3D.M[0,2]*Z

    cv.WarpAffine(src, dest, warp)

    return numpy.asarray(dest)



if __name__ == '__main__':

    directory = "/home/sphan/data/fhv6"
    basename = "fhv6a"

    p = txbr.TxBRproject( directory,basename )
    p.load()

    indexOfSeries = 0

    project(p,indexOfSeries)