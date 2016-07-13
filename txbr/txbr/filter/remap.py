import os
import numpy
import txbr
import txbr.utilities
import util
import remapFA

log = txbr.log

def remap( tracks, mask, XYZ, M0, u, angles, order, localMagnification, center= None, \
           doPlot=True, directory=None, basename=None, model="" ):
    '''Calculate the 2D remap functions to apply on the micrographs.
    ntilts is the total number of tilts in the series
    n the number of points used to calculate the remap
    Variable tracks: is a numpy array of shape (ntilts,n,2) containing 
    the locations of the tracks to strengthen.
    Variable mask:  is a numpy array of shape (ntilts,n) specifying is
    a track should be used to calculate the remap.
    Variable XYZ: is a numpy array of shape (n,3) and contains the 3D 
    estimated locations of the tracks.
    Variable M0: is a numpy array of shape (3) and is a point that belongs
    to the rotation axis of the tilt series.
    Variable u: is a numpy array of shape (3) and is a unit vector along 
    the rotation axis.
    Variable angles: is a numpy array of shape (ntilts) that specifies the
    orientation of the specimen relatively to the rotation axis (M0,u).Should
    be entered in degrees.
    Variable order: is a numpy array of shape (ntilts) that specifies the
    orientation of the specimen relatively to the rotation axis (M0,u).
    Variable localMagnification: specifies the oversampling that will be used.
    Variable center: represents the center of the image.  It will be used to
    rotate the "ideal tracks" on to the horizontal. If center is None, the 
    center of mass of XYZ will be used
    '''
    
    log.info('Calculate remapping coeffients for the micrograph stack.')
    log.info('M0: [%.1f,%.1f,%.1f]' %tuple(M0))
    log.info('u: [%.1f,%.1f,%.1f]' %tuple(u))
    
    ntilts = tracks.shape[0]
    n = tracks.shape[1]
    
    orders = numpy.repeat(order,2)
    
    if center==None:
        center = numpy.mean(XYZ,axis=0)
        
    R0 = util.Rotation3D(center,[0.0,0.0,1.0],-numpy.arccos(u[1]/numpy.dot(u,u)))
    
    data_in = []
    data_out = []
    map2D = []
    map2Dinv = []
        
    for itilt in range(ntilts):
        
        indices = numpy.where(mask[itilt,:])
        
        input = tracks[itilt,:,:]
        input = input[indices]

        XYZ_ = XYZ[indices]
        
        rotation = util.Rotation3D(M0,u,numpy.radians(angles[itilt]))
        
        base = rotation.forward(XYZ_)
        base = R0.forward(base) # Force the traces to be horizontal
        base = numpy.delete(base,2,1)
        
        # Forward Map
    
        remap_transform = util.polynomial_regression( input, base, orders )
        
        output = numpy.column_stack((remap_transform[0].eval(input),remap_transform[1].eval(input)))
        
        data_in.append(input)
        data_out.append(output) # Output as calculated by the polynomial transformation...
        
        row1 = remap_transform[0].coeff
        row2 = remap_transform[1].coeff
        row3 = numpy.zeros_like(row1)
        row4 = numpy.zeros_like(row1)
        
        row3[0] = 1.0
        row4[0] = 1.0
                
        remap = numpy.row_stack((row1,row2,row3,row4))
        
        map2D.append(remap)
        
        # Reverse Map
    
        remap_transform_inv = util.polynomial_regression( base, input, orders )
        
        output = numpy.column_stack((remap_transform_inv[0].eval(input),remap_transform_inv[1].eval(input)))
          
        row1 = remap_transform_inv[0].coeff
        row2 = remap_transform_inv[1].coeff
        row3 = numpy.zeros_like(row1)
        row4 = numpy.zeros_like(row1)
        
        row3[0] = 1.0
        row4[0] = 1.0
                
        remapinv = numpy.row_stack((row1,row2,row3,row4))
        
        map2Dinv.append(remapinv)
                
    map2D = numpy.asarray(map2D)
    map2Dinv = numpy.asarray(map2Dinv)
    
    if doPlot:
        
        txbr.utilities.plotRemap( data_in, data_out, directory, basename, model )
        
    if directory!=None and basename!=None:
        
        filename = os.path.join(directory, '%s.remap' %basename)

        log.info("Store the 2D Remap coefficients into file %s" %filename)
        
        store_inverse_map = False
        
        if store_inverse_map:
            log.warning("Stores the reverse remap!!!")
            remapFA.store2DRemap( map2Dinv, order, filename, localMagnification)
        else:
            remapFA.store2DRemap( map2D, order, filename, localMagnification)
    
    