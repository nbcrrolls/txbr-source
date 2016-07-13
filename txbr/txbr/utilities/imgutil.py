import os.path
import scipy
import numpy
import mrc
import util
import time

from PIL import Image

def array2int( u, dens = 6.0 ):
    '''Change type and normalize an array to adpat image format'''

    m1 = numpy.mean(u)
    m2 = numpy.var(u)

#    print "shape: %s" %(str(u.shape))
#    print "m1: %f" %m1
#    print "m2: %f" %m2

    u = ((u[:,::-1].T-m1)*255.0/numpy.sqrt(m2)/dens + 255.0/2.0)

    u = numpy.where(u>=255.0,255.0,u)
    u = numpy.where(u<=0.0,0.0,u)
    
    u = u.astype(numpy.uint8)

    return u


def saveAsPNG( u , name, size=(500,500) ):
    '''Save a numpy array as an png image thumbnail'''

    u = array2int(u)

    filename = "%s.gif" %name
    
    try:
        im= Image.fromarray(u,'L')
        im.thumbnail(size, Image.ANTIALIAS)
        im.save(filename, "GIF")
    except IOError:
        print "cannot create thumbnail for", filename


def make_animated_tilt_series( basename, directory, work_directory=None, durationRatio=0.05, MAX_SIZE=350, NUMBER_OF_FRAMES = 31 ):

    if work_directory==None: work_directory = directory

    mrcVolume = os.path.join( directory, "%s.preali" %basename )

    if not os.path.exists(mrcVolume):
        print "File %s does not exist!" %mrcVolume
        return

    mrcVolume = mrc.MRCFile(mrcVolume)

    nx = mrcVolume.nx
    ny = mrcVolume.ny
    nz = mrcVolume.nz

    stepz = int(nz/NUMBER_OF_FRAMES)
    bin = int(max(nx,ny)/MAX_SIZE)

    frames = []

    for iz in range(nz-1,-1,-stepz):
        print "Slice z #%i" %iz
        u = mrcVolume.getZSliceAt(iz)
        if bin!=1: u = u[::bin,::bin]
        u = array2int( u, dens = 5.0 )
        frames.insert(0,u)
        frames.append(u)
        
    file_out = os.path.join( work_directory,"%s.tilt.gif" %basename )

    util.writeGif( file_out, frames, duration=durationRatio*stepz, dither=0)



#def prepareCCDBSeries( basename, directory, work_directory ):
#    '''This routine generate an animated gif for a tilt series (in the z direction),
#    and creates a thumbnail for the middle slice'''
#
#    if work_directory==None:
#        work_directory = directory
#
#    # Make an animated gif for the preali file
#
#    make_animated_tilt_series( basename, directory, work_directory )
#
#    # Create a thumbnail from the middle slice of the preali
#
#    file_in = os.path.join( directory, "%s.st" %basename )
#    file_out = os.path.join( work_directory, "%s.thumbnail" %basename )
#
#    if os.path.exists(file_in):
#        mrcVolume = mrc.MRCFile(file_in)
#        saveAsPNG( mrcVolume.getMiddleSliceInZ(), file_out )



def make_animated_gif( filename, directory, normalizeEachSlice=False, durationRatio=0.05, 
                       showLogo=True, doXAnim=True, doYAnim=True, doZAnim=True, bin=1 ):

    doPublish = doXAnim and doYAnim and doZAnim

    t0 = time.time()
    
    mrcVolume = os.path.join( directory, filename )
    mrcVolume = mrc.MRCFile(mrcVolume)
    
#    if not normalizeEachSlice:
#        mrcVolume.updateHeader()
#        m0,m1,m2,min,max = mrcVolume.stats()

    t1 = time.time()

    print "stats: %f" %(t1-t0)
    
    nx = mrcVolume.nx
    ny = mrcVolume.ny
    nz = mrcVolume.nz

    print "(nx,ny,nz): (%i,%i,%i)" %(nx,ny,nz)


    slices = numpy.empty((nx,ny,nz))

    for iz in range(nz):
        slices[:,:,iz] = mrcVolume.getZSliceAt(iz)


    m1 = numpy.mean(slices)
    m2 = numpy.std(slices)
#    min_ = numpy.min(slices)
#    max_ = numpy.max(slices)

    dens = 10.0
    slices = ((slices-m1)*255.0/m2/dens + 255.0/2.0)

    if showLogo:
        data_dir =  os.path.join(os.path.dirname(__file__),'..','..','..','data')
        logo = os.path.join(data_dir,'img','txbr_on_the_fly_logo.png')
        image = Image.open(logo)
        arr = scipy.asarray(image)
        (nx_l,ny_l,nc) = arr.shape
        arr = (255 - arr)/2.0
        for iz in range(nz):
            slices[-nx_l:,-ny_l:,iz] += arr[:,:,0]

    slices = numpy.where(slices>=255.0,255.0,slices)
    slices = numpy.where(slices<=0.0,0.0,slices)
    slices = slices.astype('int8')


 #   slices = ((slices-min_)*255.0/(max_-min_))

    t2 = time.time()

#    print numpy.mean(slices)
#    print numpy.std(slices)

    print "load: %f" %(t2-t1)


    stepx = 4
    stepy = 4
    stepz = 2
    
    dens = 2.5

    t2 = time.time()


    if doZAnim: # Create animation in the z direction

        data_dir =  os.path.join(os.path.dirname(__file__),'..','..','..','data')
        logo = os.path.join(data_dir,'img','txbr_on_the_fly_logo.png')


    
        frames = []

        for iz in range(nz-1,-1,-stepz*bin):
#            u = mrcVolume.getZSliceAt(iz)
#            if bin!=1: u = u[::bin,::bin]
#            if normalizeEachSlice:
#                m1 = numpy.mean(u)
#                m2 = numpy.var(u)
#                print "dede"
#            u = ((u[:,::-1].T-m1)*255.0/numpy.sqrt(m2)/dens + 255.0/2.0)
#            #u = (u[:,::-1].T-min)*255.0/(max-min)
#            if showLogo: u[-nx_l:,-ny_l:] += arr[:,:,0]
##            u = numpy.where(u>=255.0,255.0,u)
##            u = numpy.where(u<=0.0,0.0,u)
#            u = u.astype(numpy.uint8)
          #  u = slices[:,:,iz].astype(numpy.uint8)
            u = slices[:,:,iz]
            frames.insert(0,u)
            frames.append(u)

        t3 = time.time()

        print "z: %f" %(t3-t2)

        util.writeGif( os.path.join( directory,'Z-anim.gif' ), frames, duration=durationRatio*stepz, dither=0)
    
    
    if doXAnim: # Create animation in the X direction


        t4 = time.time()

        frames = []

        for ix in range(nx-1,-1,-stepx):
#            u = mrcVolume.getXSliceAt(ix)
#            if bin!=1: u = u[::bin,:]
#            u = ((u[::-1,:]-m1)*255.0/numpy.sqrt(m2)/dens + 255.0/2.0)
#            #u = (u[::-1,:]-min)*255.0/(max-min)
#            u = numpy.where(u>=255.0,255.0,u)
#            u = numpy.where(u<=0.0,0.0,u)
#            u = u.astype(numpy.uint8)
           # u = slices[ix,:,:].astype(numpy.uint8)
            u = slices[ix,:,:]
            frames.insert(0,u)
            frames.append(u)

        t5 = time.time()


        print "x: %f" %(t5-t4)

        util.writeGif( os.path.join( directory,'X-anim.gif' ), frames, duration=durationRatio*stepx, dither=0)
    
    if doYAnim: # Create animation in the Y direction

        t6 = time.time()


        frames = []

        for iy in range(ny-1,-1,-stepy):
#            u = mrcVolume.getYSliceAt(iy)
#            if bin!=1: u = u[:,::bin]
#            u = ((u[:,::-1].T-m1)*255.0/numpy.sqrt(m2)/dens + 255.0/2.0)
#            #u = (u[:,::-1].T-min)*255.0/(max-min)
#            u = numpy.where(u>=255.0,255.0,u)
#            u = numpy.where(u<=0.0,0.0,u)
#            u = u.astype(numpy.uint8)
        #    u = slices[:,iy,:].astype(numpy.uint8)
            u = slices[:,iy,:]
            frames.insert(0,u)
            frames.append(u)

        t7 = time.time()


        print "x: %f" %(t7-t6)

        util.writeGif( os.path.join( directory,'Y-anim.gif' ), frames, duration=durationRatio*stepy, dither=0)

    if doPublish:

        f = open(os.path.join( directory, 'result.html' ),'wb')

        f.write('<html>')
        f.write('<head>')
        f.write('<meta http-equiv="refresh" content="10">')
        f.write('</head>')
        f.write('<body>')
        f.write('<div align="center">')
        f.write('<table align="center" cellspacing="10">')
        f.write('<tr><td> <span style="display: block;text-align: center">')
        f.write('<img src="Y-anim.gif" style="border: 0px solid black"/></span></td></tr>')
        f.write('<tr><td><span style="display: block; text-align: center">')
        f.write('<img src="Z-anim.gif" style="border: 0px solid black"/></span></td>')
        f.write('<td align="left"> <span style="display: block; text-align: center">')
        f.write('<img src="X-anim.gif" style="border: 0px solid black" /></span></td></tr></table>')
        f.write('</body>')
        f.write('</html>')

    tf = time.time()

    print "total %f" %(tf-t0)


    
