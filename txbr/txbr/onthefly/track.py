import os
import sys
import numpy
import numpy.fft
import pylab
import mrc
import modl
import phasecorr
import util
import scipy.ndimage.filters

def getParticleTemplate( size, n ):

    alpha = -1000.0
    #n = 2*size+1

    template = numpy.zeros((n,n))

    x,y = numpy.meshgrid(range(n),range(n))
    
    x = x - size
    y = y - size
    r = numpy.sqrt(x**2+y**2)

    template = numpy.where( x**2+y**2<=size**2/4.0, alpha, 0.0 )
   # template = alpha*(numpy.tanh(-500.0*(r/size-1.0/2.0)) + 1.0)


    return template


def removeBackground( u ):

    nx,ny = u.shape

#    u = scipy.ndimage.filters.gaussian_filter(u,1.2)
  #  u = scipy.ndimage.filters.maximum_filter(u,3)

    ratio = 0.20
    lx = int(ratio*nx/2.0)
    ly = int(ratio*ny/2.0)

    background = (numpy.sum(u)-numpy.sum(u[lx:-lx,ly:-ly]))/(nx*ny-(nx-2.0*lx)*(ny-2.0*ly))

    #print "(lx,ly)=(%i,%i) min,max=%.2f,%.2f  background=%.2f" %(lx,ly,numpy.min(u),numpy.max(u),background)

    th = background + 0.7*(numpy.min(u)-background)

    v = numpy.where(u<=th,u-th,0.0)

    v = scipy.ndimage.filters.gaussian_filter(v,1.2)

    return v


def centroid( u ):

    nx,ny = u.shape

    v =  removeBackground(u)

    x,y = numpy.meshgrid(range(nx),range(ny))

    x = x.T - nx/2
    y = y.T - ny/2

    N = numpy.sum(v)
    
    centroid_x = numpy.sum(v*x)/N
    centroid_y = numpy.sum(v*y)/N

    return centroid_x,centroid_y,v


def marker( size=10, n=100 ):
    '''Act as a marker (square shape of side length given by variable size) in the
    center of an array of dimension (n,n)'''

    marker = numpy.ones((n,n))

    l1 = int((n-size)/2.0)
    l2 = int((n+size)/2.0)

    marker[l1,l1:l2] = -2.0
    marker[l2,l1:l2] = -2.0
    marker[l1:l2,l1] = -2.0
    marker[l1:l2,l2] = -2.0

    return marker



def refineWithBeadModel( basename='test', model='test_test.mod', size=5, ntrack=None ):

#    command1 = "findbeads3d -size %i -thresh 0.45 test_vol.mrc beads.mod" %size
#    os.system(command1)

    stck = mrc.MRCFile("%s.st" %basename)
    nx,ny = stck.nx,stck.ny
    ntilt = stck.nz

    sx,sy = stck.sx,stck.sy
    print "scale: (%.1f,%.1f)" %(sx,sy)

    size = numpy.ceil(size*10.0/sx).astype('int')

    fiducialModel = modl.Model(model)
    fiducialModel.loadFromFile()

    mainObj = fiducialModel.objects[0]

    #ntrack = 20

    if ntrack==None:
        ntrack = len(mainObj.contours)
    else:
        ntrack = min(ntrack, len(mainObj.contours))

    mainObj.contours = mainObj.contours[:ntrack]

    points = mainObj.points()
    zmin = numpy.min(points[:,2])
    zmax = numpy.max(points[:,2])

    l = 4*size+1
    l = 6*size+1
    kernel_size = 3

    u = numpy.zeros((ntilt,l,ntrack,l)) # input
    v = numpy.zeros((ntilt,l,ntrack,l)) # output
    skip = numpy.zeros((ntilt,ntrack))
    ref = numpy.zeros(ntrack,dtype='int')
    mark = marker(2*size+1,l)

    use_centroid = True

    # Initialization
    
    for itilt in range(ntilt):
        print itilt
        slice = stck.getZSliceAt(itilt)
        for itrack in range(ntrack):
            point = mainObj.contours[itrack].getPointAtZ(itilt)
            if point.size<3:
                skip[itilt,itrack] = True
                continue
            xref1 = int(point[0,0])-l/2
            yref1 = int(point[0,1])-l/2
            xref2 = xref1+l
            yref2 = yref1+l
            skip[itilt,itrack] = xref1<0 or yref1<0 or xref2>=nx or yref2>=ny
            if not skip[itilt,itrack]:
                u[itilt,:,itrack,:] = slice[xref1:xref2,yref1:yref2]

    for itrack in range(ntrack):
        index = numpy.argwhere(skip[:,itrack]==0)
        if len(index)==0: continue
        itiltref = int(len(mainObj.contours[itrack].points)/2)
        ref[itrack] = index[numpy.abs(index-itiltref).argmin()]
        print "Track #%i: npts=%i   tiltref=%i" %(itrack,len(mainObj.contours[itrack].points),ref[itrack])

    # Tracking part

    translations = numpy.zeros((ntilt,ntrack,2))

    for itrack in range(ntrack):
        itiltref = ref[itrack]
        uref = u[itiltref,:,itrack,:]
        uref = scipy.signal.medfilt(uref,kernel_size=kernel_size)
        for itilt in range(ntilt):
            if skip[itilt,itrack]: continue
            #u1 = scipy.signal.medfilt(u[itilt,:,itrack,:],kernel_size=kernel_size)
            u1 = u[itilt,:,itrack,:]
            if use_centroid:
                tx, ty, u2 = centroid( u1)
                u[itilt,:,itrack,:] = u2
            else:
                u1 = removeBackground(u1[kernel_size:-kernel_size,kernel_size:-kernel_size])
                #u2 = getParticleTemplate(size,l)
                u2 = removeBackground(uref[kernel_size:-kernel_size,kernel_size:-kernel_size])
                tx,ty,widthx,widthy,theta = phasecorr.correlate(u1,u2)
            print "Track #%i: tilt=%i/%i  t=(%+.2f,%+.2f)  module:%.2f" %(itrack,itilt,itiltref,tx,ty,numpy.sqrt(tx**2+ty**2))
            translations[itilt,itrack,0] = tx
            translations[itilt,itrack,1] = ty


    for itilt in range(ntilt):
        slice = stck.getZSliceAt(itilt)
        for itrack in range(ntrack):
            if skip[itilt,itrack]: continue
            points = mainObj.contours[itrack].points
            if translations[itilt,itrack,0]**2 + translations[itilt,itrack,1]**2>4.0*size**2:
                print "too large"
            point = mainObj.contours[itrack].getPointAtZ(itilt)

            point[0,0] = int(point[0,0]) + translations[itilt,itrack,0]
            point[0,1] = int(point[0,1]) + translations[itilt,itrack,1]
            point = numpy.squeeze(point)

            index = numpy.argwhere(points[:,2]==itilt).ravel()
            for ipt in index:
                mainObj.contours[itrack].points[ipt,0] = point[0]
                mainObj.contours[itrack].points[ipt,1] = point[1]

            if point.size!=3: continue
            print "Track #%i: tilt=%i  pt=(%+.2f,%+.2f)" %(itrack,itilt,point[0],point[1])
            xref1 = int(point[0])-l/2
            yref1 = int(point[1])-l/2
            xref2 = xref1+l
            yref2 = yref1+l
            if xref1<0 or yref1<0 or xref2>=nx or yref2>=ny: continue
            v[itilt,:,itrack,:] = removeBackground(slice[xref1:xref2,yref1:yref2])*mark

    fiducialModel.filename = "%s-otf.fid" %basename

    fiducialModel.save()

    output = mrc.MRCFile("%s-otf.mrc" %basename)
    output.setHeader(ntilt*l,ntrack*l,2)
    output.setZSliceAt(0,numpy.resize(u,(ntilt*l,ntrack*l)))
    output.setZSliceAt(1,numpy.resize(v,(ntilt*l,ntrack*l)))
    output.updateHeader()

    os.system("imod %s-otf.mrc" %basename)

    


if __name__=='__main__':


    refineWithBeadModel(basename='FHV-16_gatan', model='FHV-16_gatan_test.mod', size=5, ntrack=None)
    sys.exit()



    command = "findbeads3d -size 5 -thresh 0.45 test_vol.mrc beads.mod"

    filename = '/ncmirdata3/sphan/FHV-14/de12/bin4/de12.st'
    f = mrc.MRCFile(filename)

    print f.sx
    print f.sy

    fout = mrc.MRCFile("test.out",template=filename)

    nz = f.nz

    size = int(50.0/f.sy) + 4.0
    
    template = getParticleTemplate( size )

    pylab.matshow(template)
    pylab.show()

    for iz in range(nz):

        u = f.getZSliceAt(iz)

        v = util.crossCorrelate( u, template, crop_rate=0.0, ret='array' )
    
        th = numpy.mean(v) + 0.3*(numpy.max(v)-numpy.mean(v))
        #v = numpy.where(v>th,100.0,0.0)
        v = numpy.where(v>th,1000.0*(v+th),0.0)

        res = numpy.zeros_like(u)
        res[size:-size,size:-size] = v

#        array2d_m = scipy.ndimage.filters.maximum_filter(res,3)
#        res = ((res==array2d_m) & (res!=0))
#        v = 1000.0*v
        #v = numpy.where(v>th,v,0.0)

        n = len(numpy.argwhere(res!=0.0))

        print "iz=%i  n=%i" %(iz,n)

        fout.setZSliceAt(iz,res)

    fout.updateHeader()

#    pylab.matshow(v)
#    pylab.show()

#    mrc.createMRCFile("test.out",res,True)