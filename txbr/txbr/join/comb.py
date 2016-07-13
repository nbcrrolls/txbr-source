import os
import numpy
import mrc
import txbr
import re
import util

from txbr import log
from txbr import X_TAPERING_RATIO
from txbr import Y_TAPERING_RATIO
from txbr.txbrdao import loadProject

# Layering options
MASK_MERGING = 0
MASK_LAYERING = 1


def getBlockingMaskFromFile( basename, bin=None ):

    if bin==None: bin=1.0
    else: bin = float(bin)

    mask = {}

    filename = "%s.%s" %(basename,"mask")

    if os.path.exists(filename):

        f = open(filename)

        for line in f:

            m = re.match('(\S*)->\((\S*):(\S*),(\S*):(\S*)\)', line)
            if not m: continue

            extension = m.group(1)
            xstart,xstop = float(m.group(2)),float(m.group(3))
            ystart,ystop = float(m.group(4)),float(m.group(5))

            mask[basename + extension] = [[xstart/bin,xstop/bin],[ystart/bin,ystop/bin]]

        f.close()

    return mask


def writeBlockingMaskToFile( basename, data, bin=None ):

    if bin==None: bin=1.0
    else: bin = float(bin)
    
    filename = "%s.%s" %( basename, "mask" )

    f = open(filename,"w")

    for extension,ll,ur in data:
        f.write("%s->(%i:%i,%i:%i)\n" %(extension,ll[0]*bin,ur[0]*bin,ll[1]*bin,ur[1]*bin))

    f.close()

    return



def getFinalFrame( basenames, LLcorners, URcorners, bin=None ):
    '''Returns the final frame.
    Variable LLcorners: a numpy array of shape (n,2) containing
    the lower left points for each mosaic tile reconstruction.
    Variable URcorners: a numpy array of shape (n,2) containing
    the upper right points for each mosaic tile reconstruction.
    '''
    
    # Get the final frame width and height
    
    bottom_left = numpy.min(LLcorners,axis=0)
    top_right = numpy.max(URcorners,axis=0)
    
    xstart,ystart,zstart = bottom_left
    
    nx = top_right[0] - bottom_left[0] + 1
    ny = top_right[1] - bottom_left[1] + 1
    nz = top_right[2] - bottom_left[2] + 1
    
    # Size and offset of each mosaic
    
    offsets = LLcorners - bottom_left
    sizes = URcorners - LLcorners + 1

    # Get the additional mask offsets and size applied during the stitching

    LLcorners_mask = LLcorners.copy()
    URcorners_mask = URcorners.copy()

    basename, extensions = txbr.utilities.extract_series_name( basenames )
    blkmask = getBlockingMaskFromFile( basename, bin=bin )

    for index,name in enumerate(basenames):
        if name in blkmask.keys():
            x,y = blkmask[name]
            LLcorners_mask[index,0] = x[0]
            LLcorners_mask[index,1] = y[0]
            URcorners_mask[index,0] = x[1]
            URcorners_mask[index,1] = y[1]

    offsets_mask = LLcorners_mask - bottom_left
    sizes_mask = URcorners_mask - LLcorners_mask + 1
    
    # Some display
        
    for index,(name,offset,size) in enumerate(zip(basenames,offsets,sizes)):
        log.info("Mosaic #%i (%s): offset: %s  size: %s" %(index,basenames[index],offset,size))

    writeBlockingMaskToFile( basename, sorted(zip(extensions,LLcorners_mask,URcorners_mask)), bin=bin )

    return (xstart,ystart,zstart), (nx,ny,nz), offsets, sizes, offsets_mask, sizes_mask


def trimMasks( masks, delta_xl, delta_xr, delta_yb, delta_yt ):

    for index,m in enumerate(masks):
        m[:delta_xl[index],:] = 0.0
        m[-delta_xr[index]:,:] = 0.0
        m[:,:delta_yb[index]] = 0.0
        m[:,-delta_yt[index]:] = 0.0
        
        
def smoothMasks( masks, delta_xl, delta_xr, delta_yb, delta_yt, lx, ly ):

    for index,m in enumerate(masks):
        for wx in range(lx[index]):
            m[delta_xl[index]+wx,:] *= (wx+1.0)/lx[index]
            m[-delta_xr[index]-wx,:] *= (wx+1.0)/lx[index]
        for wy in range(ly[index]):
            m[:,delta_yb[index]+wy] *= (wy+1.0)/ly[index]
            m[:,-delta_yt[index]-wy] *= (wy+1.0)/ly[index]
        

def getMask( (nx,ny), offsets, sizes, mode=MASK_MERGING, offsets_blkmask=None, sizes_blkmask=None, bin=None ):
    '''Calculate a mask.
    Variable offsets: Array containing the postion (in 2D) of the bottom left corner for each tile
    Variable sizes: Array containing the size (in 2D) of each tile
    Variable offsets_blkmask: Similar to variable offsets for the all or nothing mask
    Variable sizes_blkmask: Similar to variable sizes for the all or nothing mask
    '''
        
    mask = numpy.zeros((nx,ny)) # full mask
    masks = [ numpy.ones((size[0],size[1])) for size in sizes ]
    
    delta_x = sizes[:,0]*X_TAPERING_RATIO
    delta_y = sizes[:,1]*Y_TAPERING_RATIO
    
    delta_xl = delta_x
    delta_xr = delta_x
    delta_yb = delta_y
    delta_yt = delta_y

    delta_LL = offsets_blkmask - offsets
    delta_UR = offsets + sizes - offsets_blkmask - sizes_blkmask

    delta_xl = numpy.where(delta_LL[:,0]>=delta_x,delta_LL[:,0],delta_x)
    delta_xr = numpy.where(delta_UR[:,0]>=delta_x,delta_UR[:,0],delta_x)
    delta_yb = numpy.where(delta_LL[:,1]>=delta_y,delta_LL[:,1],delta_y)
    delta_yt = numpy.where(delta_UR[:,1]>=delta_y,delta_UR[:,1],delta_y)

    trimMasks( masks, delta_xl, delta_xr, delta_yb, delta_yt )
    
    if mode==MASK_MERGING: # Same weight for tiles in common part
        
        lx,ly = 2*delta_x,2*delta_y
        
        lx = numpy.ceil(lx).astype('int')
        ly = numpy.ceil(ly).astype('int')
        
        smoothMasks( masks, delta_xl, delta_xr, delta_yb, delta_yt, lx, ly )
    
        for index,(offset,size) in enumerate(zip(offsets,sizes)):
            mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]] += masks[index][:,:]
                    
        for index,(offset,size) in enumerate(zip(offsets,sizes)):
            masks[index][:,:] = masks[index][:,:]/mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]]
        
    elif mode==MASK_LAYERING:
        
        for index,(offset,size) in enumerate(zip(offsets,sizes)):
            mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]] += masks[index][:,:]
            masks[index][:,:] = mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]]
            masks[index][:,:] = numpy.where(masks[index][:,:]==2,0,masks[index][:,:])
            mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]] = numpy.where(mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]]==2,1,mask[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]])
                
    trimMasks( masks, delta_xl, delta_xr, delta_yb, delta_yt )
                
    return masks


def combineVolumes( mrc_inputs, mrc_output, start_output, size_output, masks, offsets, sizes, slices=None):

    n = offsets.shape[0]
    
    (xstart,ystart,zstart) = start_output
    (nx,ny,nz) = size_output
    
    if slices==None:
        nz = mrc.MRCFile(mrc_inputs[0]).nz
        slices = range(nz)
    else:
        nz = len(slices)
        
    final_vol = mrc.MRCFile(mrc_output)
    final_vol.setHeader(nx,ny,nz)
    
    mrc_inputs = [ mrc.MRCFile(file) for file in mrc_inputs ]
    
    array = numpy.zeros((nx,ny))
  
    for index,iz in enumerate(slices):

        log.info("joining slice #%i" %iz)
        
        array[:,:] = 0
        
        for i in range(n):
            offset = offsets[i]
            size = sizes[i]
            slice = mrc_inputs[i].getZSliceAt(iz)
            array[offset[0]:offset[0]+size[0],offset[1]:offset[1]+size[1]] += masks[i]*slice
        
        final_vol.setZSliceAt(index,array)
        
    final_vol.updateGrid( nxstart=xstart, nystart=ystart, nzstart=zstart, mx=nx, my=ny, mz=nz )
    final_vol.updateHeader()
    
    
def joinVolumes( directory, basenames, work_directory, output=None, check_mask=False, zmin=None, zmax=None, test=False, bin=None ):
    '''Join Series Volumes together.'''
    
    log.info("Join Series Volumes!")
     
    if isinstance( basenames, str ):
        basenames = txbr.utilities.extract_series_name_from_dir( basenames, work_directory, extension=txbr.utilities.TXBR_EXTENSION )
                
    LLcorners = []
    URcorners = []
    mrc_inputs = []
    
    for basename in basenames:

        project = loadProject(directory=directory,basenames=basename,work_directory=work_directory, bin=bin)
       
        origin = project.reconstruction.origin
        end = project.reconstruction.end

        if zmin!=None and zmax!=None:
            LLcorners.append((origin.x,origin.y,zmin))
            URcorners.append((end.x,end.y,zmax))
            mrc_inputs.append("%s_z_%.1f.out" %(os.path.join(project.backprojection_directory,basename),zmin))
        else:
            LLcorners.append((origin.x,origin.y,origin.z))
            URcorners.append((end.x,end.y,end.z))
            mrc_inputs.append("%s_z_%.1f.out" %(os.path.join(project.backprojection_directory,basename),origin.z))

    LLcorners = numpy.asarray( LLcorners, dtype='int' )
    URcorners = numpy.asarray( URcorners, dtype='int' )
        
    LLcorners = LLcorners[:,:]
    URcorners = URcorners[:,:]
    
    log.info("Lower Left Corners: %s" %(",".join(str(LLcorners).splitlines())))
    log.info("Upper Right Corners: %s" %(",".join(str(URcorners).splitlines())))

    (xstart,ystart,zstart), (nx,ny,nz), offsets, sizes, offsets_blkmask, sizes_blkmask = getFinalFrame( basenames, LLcorners, URcorners, bin=bin )

    masks = getMask( (nx,ny), offsets, sizes, offsets_blkmask=offsets_blkmask, sizes_blkmask=sizes_blkmask, bin=bin )
    
    if output==None: # build a standart name for the output
        basename, extensions = txbr.utilities.extract_series_name(basenames)
        output = "%s_z_%.1f.out" %(basename ,zstart)
        
    if check_mask:
        msk_file = mrc.MRCFile(output)
        msk_file.setHeader(nx,ny,1)
        msk_file.setZSliceAt(0,masks)
        msk_file.updateHeader()
        os.system("imod %s" %msk_file)
        return
    
    if test:
        slices = [ int(nz/2) ]
    else:
        slices = None
        
    combineVolumes( mrc_inputs, output, (xstart,ystart,zstart), (nx,ny,nz), masks, offsets, sizes, slices )
    
    if test:
        os.system('imod %s' %output)
        
        
def weightVolume( volume_file, weigth_file, output, test=False ):
    '''
    Variable volume_file: The name of the MRC file to weight
    Variable weigth_file: The name of the MRC file that specify the weight of 
    each pixel in the volume. The weight file can be generated by backprojecting
    slices made of 1.
    '''
    
    vol = mrc.MRCFile( volume_file )
    target = mrc.MRCFile( weigth_file )
    
    if vol.nx!=target.nx or vol.ny!=target.ny or vol.nz!=target.nz:
        return;
    
    nx, ny, nz = vol.nx, vol.ny, vol.nz
    
    slices = range(nz)
    if test:
        slices = [ int(vol.nz/2) ]
    
    dest = mrc.MRCFile( output )
    dest.setHeader( nx, ny, nz )
    
    dest_slice = numpy.zeros((nx,ny))
    
    for index,iz in enumerate(slices):
        slice = vol.getZSliceAt(iz)
        tg = target.getZSliceAt(iz)
        dest_slice = numpy.where(tg!=0,slice/tg,slice)
        dest.setZSliceAt( index, dest_slice )
        
    dest.updateHeader()
        
        
def joinMicrographs( directory, basenames, work_directory, output=None, test=True ):
    '''Join micrographs together. In this routine, the 3D alignment of the mosaics is used to generate a
    stack of 2D stitched micrographs. This is useful in the case of the 4000#2 microscope with its four cameras.'''
    
    log.info("Join micrographs together!")
     
    if isinstance( basenames, str ):
        basenames = txbr.utilities.extract_series_name_from_dir( basenames, work_directory, extension=txbr.utilities.TXBR_EXTENSION )

    series = []    # The set of series; each series is loaded individually
    rot_axes = []    # Rotation Axis for each series
    projections = []    # Projection map for every series
    dest_frames = []    # Destination Frame for the zero-tilt slice of each series
        
    for basename in basenames:
        project = loadProject(directory=directory,basenames=basename,work_directory=work_directory)
        series = series + project.series
        
    if len(series)==0: 
        og.error("No series available!")
        return;
    
    angles = series[0].tiltAngles
        
    for index,s in enumerate(series):
        log.info("Series Basename: %s" %s.basename)
        log.info("(nx,ny): (%i,%i)" %s.dimensions())
        iref0 = s.indexOfReferenceExposure()
        M = [ s.projection.x_coefficients[iref0], \
              s.projection.y_coefficients[iref0], \
              [0.0,0.0,0.0,1.0] ]
        M = numpy.row_stack(M)
        # Include the preali transformation ? to check
        M[0,0] -= s.prealignTranslations[iref0,0]
        M[1,0] -= s.prealignTranslations[iref0,1]
        tfinv = util.AffineTransformation3D(M).inv()
        src_frame = numpy.row_stack((s.getSourceFrame(),[0.0,0.0,0.0,0.0])).T
        dest_frame = tfinv.forward(src_frame)
        dest_frames.append(dest_frame)
        rot_axes.append(s.rotAxis)
        
    dest_frames = numpy.row_stack(dest_frames)
        
    bottom_left = numpy.min(dest_frames,axis=0)
    top_right = numpy.max(dest_frames,axis=0)
    
    corners = numpy.row_stack(( bottom_left, top_right )).astype('int')
    middle = numpy.average( corners, axis=0 )    # Rotation axis point
    
    nx = corners[1,0] - corners[0,0] + 1
    ny = corners[1,1] - corners[0,1] + 1
    
    rot_axes = numpy.average( rot_axes, axis=0 )
    
#    print "Final stack: [(nx,ny)=(%i,%i)]" %( nx, ny )
#    print "Corners: LL [%.2f,%.2f]     UR [%.2f,%.2f]" %( corners[0,0], corners[0,1], corners[1,0], corners[1,1] )
#    print "Rot Axis: Point (%.2f,%.2f,%.2f)    Direction (%.2f,%.2f,%.2f) " %( middle[0], middle[1], middle[2], rot_axes[0], rot_axes[1], rot_axes[2] )
    
    nz = len(angles)
        
    if test:
        slices = range(iref0,iref0+2)
    else:
        slices = range(nz)
    
    f = mrc.MRCFile("test.mrc")
    f.setHeader( nx, ny, len(slices) )    
    
    BORDER = 0.01*nx    # Don't take the border
    
    for islice,slice in enumerate(slices):
        
        angle = angles[slice]

#        print "Index: %3i      Slice #%3i     Angle: %5.2f" %( islice, slice, angle )
        
        rot = util.Rotation3D(middle, rot_axes, -angle, in_radians=False)
        
        Gt = rot.T
        Gm = rot.M
            
        G = numpy.column_stack((Gt,Gm))
        #G[2,0] = 0.0
        #G[2,1] = 0.0
        #G[2,2] = 0.0
        #G[2,3] = 1.0
            
        rot = util.AffineTransformation3D(G)
        rot = util.AffineTransformation3D(numpy.eye(3,4,1))
        
        masks = []
        offsets = []
        tiles = []
        
        for index,s in enumerate(series):
            
            #print "Series Basename: %s   Angle: %f" %(s.basename,angle)
            #print "Index of reference %i" %(s.indexOfReferenceExposure())
            #print "Rotation: %s" %(util.AffineTransformation3D.__repr__(rot))
           
            f_ = mrc.MRCFile("%s.st" %basenames[index])
            u_ = f_.getZSliceAt(slice)
            M = [ s.projection.x_coefficients[slice], \
                  s.projection.y_coefficients[slice], \
                  [0.0,0.0,0.0,1.0] ]
            M = numpy.row_stack(M)
            # Make sure it is the projection map for the st stack
            M[0,0] -= s.prealignTranslations[slice,0]
            M[1,0] -= s.prealignTranslations[slice,1]
            
            tf = util.AffineTransformation3D(M)
            tfinv = tf.inv()
            
            g = rot.compose(tfinv)
                        
            T = g.T[:2].astype('int')
            
            M = numpy.column_stack((g.T[:2]-T,g.M[:2,:2]))
            ff = util.AffineTransformation2D(M)
            
            u_ = util.WarpAffine(u_,ff)
            
            mask = numpy.zeros_like(u_)
            mask[BORDER:-BORDER,BORDER:-BORDER] = 1
            
            mask = util.WarpAffine(mask,ff)            
            
            offsets.append((T[0],T[1]))
            masks.append(mask)
            tiles.append(u_)
            
        final_mask = numpy.zeros((nx,ny))
        final_slice = numpy.zeros((nx,ny))
        
        for offset,mask in zip(offsets,masks):
            
            ox_dest = max( offset[0] - corners[0,0], 0 )
            oy_dest = max( offset[1] - corners[0,1], 0 )
            
            ox_src = max( -ox_dest, 0 )
            oy_src = max( -oy_dest, 0 )
            
            nx_ = min( mask.shape[0]-ox_src, nx-ox_dest )
            ny_ = min( mask.shape[1]-oy_src, ny-oy_dest )
            
            final_mask[ox_dest:ox_dest+nx_,oy_dest:oy_dest+ny_] += mask[ox_src:ox_src+nx_,oy_src:oy_src+ny_]
            
        for offset,mask,tile in zip(offsets,masks,tiles):
            
            ox_dest = max( offset[0] - corners[0,0], 0 )
            oy_dest = max( offset[1] - corners[0,1], 0 )
            
            ox_src = max( -ox_dest, 0 )
            oy_src = max( -oy_dest, 0 )
            
            nx_ = min( mask.shape[0]-ox_src, nx-ox_dest )
            ny_ = min( mask.shape[1]-oy_src, ny-oy_dest )
            
            # BLEND
            #final_slice[ox:ox+nx_,oy:oy+ny_] += mask[:nx_,:ny_]*tile[:nx_,:ny_]/final_mask[ox:ox+nx_,oy:oy+ny_]
            # LAYER
            final_slice[ox_dest:ox_dest+nx_,oy_dest:oy_dest+ny_] = numpy.where(mask[ox_src:ox_src+nx_,oy_src:oy_src+ny_]!=0,tile[ox_src:ox_src+nx_,oy_src:oy_src+ny_],final_slice[ox_dest:ox_dest+nx_,oy_dest:oy_dest+ny_])
    
        f.setZSliceAt( islice, final_slice )
        
    f.updateHeader()
    
    os.system('imod test.mrc')
    
