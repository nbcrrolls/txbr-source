import os
import numpy
import mrc
import mrcfile

def createMRCFile( filename, u, scale=None, view=False ):
    '''Quick and dirty utility to dump a 3D numpy array into an MRC file.

    :param filename: The name of the MRC file to create.
    :param u: The numpy array containing the volume data to store.
    :param scale: The scale of this MRC file (pixel size)
    :param view: Open the volume with the 3dmod viewer.
    '''

    u = numpy.asarray(u)

    if len(u.shape)==2:
        nx,ny = u.shape
        nz = 1
        u = numpy.resize(u,(nx,ny,nz))
    else:
        nx,ny,nz = u.shape

    output = mrc.MRCFile(filename)
    output.setHeader(nx,ny,nz)
    if scale!=None: # Does not do anything
        output.sx = scale[0]
        output.sy = scale[1]
        output.sz = scale[2]
    for iz in range(nz):
        output.setZSliceAt(iz,u[:,:,iz])
    output.updateHeader()

    if view:
        os.system("imod %s" %filename)


def newstack( files_in, file_out, secs=None, bin=1, bounds3D=None, descriptions=None, clean=False):
    """Helper routine to stack several MRC files into a single MRC file.

    :param files_in: A list containing the paths of the MRC files to stack.
    :param file_out: The path of the output MRC file.
    :param secs: A list containing the slice indices to consider for each inputs.
        The list should have the same length as *files_in*. If None, all slices will
        be taken.
    :param bin: The bin factor.
    :param bounds3D: The bounding information.
    :param clean: If True the input files are removed from the file system.
    """

    files_in = [ f for f in files_in if os.path.lexists(f) ]
    mrc_in = [ mrc.MRCFile(fin) for fin in files_in ]

    if secs==None:
        secs = [ range(fin.nz) for fin in mrc_in ]

    dims = numpy.array([[fin.nx,fin.ny,fin.nz] for fin in mrc_in],dtype='int')
    
    nx = numpy.mean(dims[:,0]/bin).astype('int')
    ny = numpy.mean(dims[:,1]/bin).astype('int')
    nz = numpy.sum([len(s) for s in secs]).astype('int')

    print "nx={0:}  nx={1:}  nx={2:}".format(nx,ny,nz)

    if os.path.exists(file_out):
        os.remove(file_out)

    fout = mrc.MRCFile(file_out)
    fout.setHeader(nx,ny,nz)

    index = 0
    for fin,slices in zip(mrc_in,secs):
        for iz in slices:
            print "[%s] %s %i" %(file_out,fin.filename,iz)
            u = fin.getZSliceAt(iz)
#            print "BEGIN"
#            print u.shape
            u = u.reshape(nx, bin, ny, bin)
            u = numpy.squeeze(u.sum(axis=3).sum(axis=1))
#            print u.shape
#            print "END"
            fout.setZSliceAt(index, u)
            index += 1

    fout.updateHeader()

    if bounds3D!=None:
        xstart,ystart,zstart = int(bounds3D.origin.x),int(bounds3D.origin.y),int(bounds3D.origin.z)
        xstop,ystop,zstop = int(bounds3D.end.x),int(bounds3D.end.y),int(bounds3D.end.z)
        xscale,yscale,zscale = bounds3D.scale.x,bounds3D.scale.y,bounds3D.scale.z
        mrcfile.update_grid( file_out, xstart, ystart, zstart, xstop-xstart+1, ystop-ystart+1, zstop-zstart+1)
        mrcfile.update_scale( file_out, xscale, yscale, zscale)
        mrcfile.update_origin( file_out, (1-xstart)*xscale, (1-ystart)*yscale, -zstart*zscale)

    if descriptions!=None:
        for label in descriptions: mrc.add_label( output, label)

    if clean:
        for fin in files_in: os.remove(fin)


def translateMRCFile( input, output, translations, mode=None, Lx=500, Ly=500 ):
    '''Restack a MRC file given a set of translations.

    Intensity is adjusted between sequential slices. For this purpose, only part
    of an image is used (a centered window frame of half-size Lx and Ly in the x and y
    directions)

    :param input: The name of the original MRC file
    :param output: The name of the output MRC file
    :param translations: A numpy array of shape (n,3). The first column, contains the
        indices of the images to keep, while the second and third column are the actual
        translation vector values to correct for.
    :param mode: The mode of the new MRC File.
    :param Lx, Ly: Half-size value of the window frame used for intensity matching.
    '''

    if not os.path.exists(input): return
    if len(translations)==0: return

    txmin = numpy.floor(numpy.min(translations[:,1])).astype('int')
    txmax = numpy.ceil(numpy.max(translations[:,1])).astype('int')
    tymin = numpy.floor(numpy.min(translations[:,2])).astype('int')
    tymax = numpy.ceil(numpy.max(translations[:,2])).astype('int')

    txmin = numpy.min(txmin,0)  # To avoid problems with array indices (this case should not happen tough)
    tymin = numpy.min(tymin,0)

    print "txmin=%i  txmax=%i" %(txmin,txmax)
    print "tymin=%i  tymax=%i" %(tymin,tymax)

    fin = mrc.MRCFile(input)

    nx = fin.nx + txmax - txmin
    ny = fin.ny + tymax - tymin
    nz = len(translations)

    x1 = - txmin
    x2 = - txmin + fin.nx
    y1 = - tymin
    y2 = - tymin + fin.ny

#    print "x1:x2=%i:%i" %(x1,x2)
#    print "y1:y2=%i:%i" %(y1,y2)

    # Frame to adjust the intensity

    xx1 = max(0,fin.nx/2-Lx)
    xx2 = min(fin.nx,fin.nx/2+Lx)
    yy1 = max(0,fin.ny/2-Ly)
    yy2 = min(fin.ny,fin.ny/2+Ly)

    # print "x: %i:%i  y: %i:%i" %(x1,x2,y1,y2)

#    fout = mrc.MRCFile( output, template=input )
    fout = mrc.MRCFile( output )
    # TOOOOOOOOOOOOO-----------
#    fout.setHeader( nx, ny, nz, mode=fin.mode )

    for index,slice in enumerate(translations):
        if index <800:
            print "WARNING to remove constraint index<800"
            continue
        iz = int(round(slice[0]))
        tx,ty = slice[1:].astype('float')
        u = numpy.zeros((nx,ny))
        v = fin.getZSliceAt(iz)
        m = numpy.mean(v[xx1:xx2,yy1:yy2])
        std = numpy.std(v[xx1:xx2,yy1:yy2])
        v = (v - m)*100.0/std
        print "%i/%i: slice #%i  t=(%.1f,%.1f)  mean=%.1f  std=%.1f" %(index,nz,iz,tx,ty,m,std)
#        # ---Use the ndimage shift
#        u[x1:x2,y1:y2] = v
#        u = scipy.ndimage.shift(u,(tx,ty),order=0)  # it seems there is a bug for very large images
        #----Do a linear interpolation----
#        tx0 = numpy.floor(tx)
#        px = tx - tx0
#        ty0 = numpy.floor(ty)
#        py = ty - ty0
        tx0 = numpy.round(tx)
        ty0 = numpy.round(ty)
        x1_ = x1 + tx0
        y1_ = y1 + ty0
        x2_ = x1_ + fin.nx
        y2_ = y1_ + fin.ny
        u[x1_:x2_,y1_:y2_] += v
        #---------------------------------

        fout.setZSliceAt(index,u)

    fout.updateHeader()

    return
