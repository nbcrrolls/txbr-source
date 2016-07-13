import logging,logging.config
import scipy,numpy,numpy.linalg
import util

log = logging.getLogger()

def trace_angle(tracks,mask,axis):
    '''This routine evaluates trace angles.
    Variable tracks is a numpy array (of shape (numtlts,numtracks,2)) that
    contain track informations of point like markers in a micrograph series.
    Variable mask (a numpy array of shape (numtlts,numtracks)) allows to
    eventually remove points during the calculation.'''

    (numtlts,numtracks) = tracks.shape[:2];

    numtracks = numpy.sum(numpy.any(mask[:,:],axis=0))  # There might be some empty tracks
    
    #print "Number of tracks: %i" %numtracks

    ref_angle = numpy.arctan2(axis[1],axis[0]) - numpy.pi/2.0

    a,b = numpy.zeros(numtracks),numpy.zeros(numtracks)
    
    indices = numpy.argwhere(numpy.any(mask[:,:],axis=0)).ravel()
    
    index = 0

    #for itrack in range(numtracks):
    for itrack in indices:

        imask = numpy.argwhere(mask[:,itrack]==1)
        
#        if len(imask)==0:
#            continue

        x = numpy.squeeze(tracks[imask,itrack,0])
        y = numpy.squeeze(tracks[imask,itrack,1])
        
        # ADDED for mosaic and empty particles?
        if len(x.shape)==0:
            x = numpy.array([x])
            y = numpy.array([y])
        #
        
        a[index],b[index] = scipy.polyfit(x,y,1)
        
        index += 1

    log.debug('a = %s ...' %str(a)[:70])
    log.debug('b = %s ...' %str(b)[:70])

    slope = numpy.row_stack((numpy.arctan(a),numpy.arctan(1/a)))

    theta = slope.mean(axis=1)
    std = slope.var(axis=1)

    log.debug('theta = %s    std = %s' %(theta,std))

    if std[0]<=std[1]:
        angle = theta[0] - ref_angle
    else:
        angle = - theta[1] - ref_angle
        angle = - numpy.pi/2 + angle%numpy.pi

    if numtracks==0: angle = 0

    angle = - numpy.pi/2.0 + (numpy.pi/2.0+angle)%numpy.pi   # Set the trace angle in [-90,90]

    log.debug('ref_angle = %f    angle = %f' %(ref_angle,angle))

    return (ref_angle,angle)


def relative_rot_transform( tracks, mask, index_ref_series=0 ):
    '''This routine evaluates the best orthogonal transformation between
    a sets of tracks (aka two reference images between tilt series).
    Variable tracks: is a numpy array. Number of rows in variable tracks
    represent the number of considered sets.
    Variable mask: allows to remove points in the calculation.
    Variable index_ref_series is the index of the reference series from which transformations
    are evaluated.'''
    
    log.info("Calculate the relative affine transformations between series!")
    
    (nseries,numtracks) = tracks.shape[:2]
    
    path = numpy.zeros((nseries,nseries))
    transformations = []
    
    for iseries_ref in range(nseries):
        
        tf_ = []
        
        x_ref = tracks[iseries_ref,:]
        
        for iseries in range(nseries):
            
            imask = numpy.argwhere(mask[iseries_ref,:]*mask[iseries,:]==1)
        
            if len(imask)==0:
                #print numpy.argwhere(mask[iseries_ref,:]==1)
                #print numpy.argwhere(mask[iseries,:]==1)
                log.warning("No path between series #%i and #%i on the reference exposure!" %(iseries_ref,iseries))
                path[iseries_ref,iseries] = numpy.Inf
                tf_.append((None,None))
                continue
            
            x_ref_ = numpy.squeeze(x_ref[imask])
            x_ = numpy.squeeze(tracks[iseries,imask])
            
            ortreg = util.orthogonal_regression_2d(x_ref_,x_)
            
            theta = numpy.array([0.0,0.0,ortreg[1]])
            
            coeffs = numpy.eye(3,4,1)
            coeffs[:2,:3] = numpy.resize(ortreg[0].coeff,(2,3))
            
            path[iseries_ref,iseries] = 1
            transformation = util.AffineTransformation3D(coeffs)
            
            tf_.append((transformation,theta))
            
        transformations.append(tf_)
        
    # Find the most appropriate path between series index_ref_series and all 
    # the others. We use a Floyd Marhal algorith type
    
    next = numpy.zeros((nseries,nseries),dtype='int')
    next[:,:] = -1  # No intermediate
    
    for k in range(nseries):
        for i in range(nseries):
            for j in range(nseries):
                if path[i,k] + path[k,j] < path[i,j]:
                    path[i,j] = path[i,k]+path[k,j];
                    next[i,j] = k
 
    def GetPath (i,j):
        
        if path[i,j]==numpy.Inf:
            return None
        x = next[i,j]
        if x==-1:
            return transformations[i][j][0]
        else:
            tix = GetPath(i,x)
            txj = GetPath(x,j)
            if tix==None or txj==None:
                return None
            else:
                return txj.compose(tix)
        
    for i in range(nseries):
        for j in range(nseries):
            log.debug("Relative Transform between Series (%i,%i): %s" %(i,j,GetPath (i,j)))
    
    # Return results with index_ref_series as a reference
    
    tf_ = [ GetPath(index_ref_series,j) for j in range(nseries) ]
    ang_ = [item[1] for item in transformations[index_ref_series]]

    return ( tf_, ang_ )

