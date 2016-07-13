import numpy
import cv


def showArray( u ):
    '''Display a 2D array in a window'''
    
    cv.NamedWindow("image",cv.CV_WINDOW_AUTOSIZE)
    cv.ShowImage("image",cv.fromarray(u))
    cv.WaitKey(0)


def __loadImage( nameOfImage ):

    image =  cv.LoadImageM(nameOfImage)

    u = numpy.asarray(image)

    # Return only one channel (the image should be grayscale)

    if len(u.shape)>2:
        return u[::-1,:,0].T.astype('float')
    else:
        return u[::-1,:].T.astype('float')


def __saveImage( nameOfImage, u):

    print "Save image %s" %nameOfImage

    array = numpy.asarray(u[:,::-1].T.copy()).astype('float32')

    cv.SaveImage( nameOfImage, cv.fromarray(array) )
    

def cv2array(im):
    
    depth2dtype = { cv.IPL_DEPTH_8U: 'uint8',
                    cv.IPL_DEPTH_8S: 'int8',
                    cv.IPL_DEPTH_16U: 'uint16',
                    cv.IPL_DEPTH_16S: 'int16',
                    cv.IPL_DEPTH_32S: 'int32',
                    cv.IPL_DEPTH_32F: 'float32',
                    cv.IPL_DEPTH_64F: 'float64',
                    }
  
    arrdtype = im.depth
    a = numpy.fromstring( im.tostring(), dtype=depth2dtype[im.depth], count=im.width*im.height*im.nChannels)
    a.shape = (im.height,im.width,im.nChannels)

    return a
    
def array2cv(a):
    
    dtype2depth = { 'uint8':   cv.IPL_DEPTH_8U,
                    'int8':    cv.IPL_DEPTH_8S,
                    'uint16':  cv.IPL_DEPTH_16U,
                    'int16':   cv.IPL_DEPTH_16S,
                    'int32':   cv.IPL_DEPTH_32S,
                    'float32': cv.IPL_DEPTH_32F,
                    'float64': cv.IPL_DEPTH_64F,
                    }

    try:
        nChannels = a.shape[2]
    except:
        nChannels = 1
        
    cv_im = cv.CreateImageHeader( (a.shape[1],a.shape[0]), dtype2depth[str(a.dtype)], nChannels)
    cv.SetData(cv_im, a.tostring(), a.dtype.itemsize*nChannels*a.shape[1])
    
    return cv_im


def crossCorrelate( array1, array2, crop_rate=0.1, ret='max' ):
    '''Cross-correlate two numpy array by changing them into openCV MAT objects'''

    array1 = numpy.asarray(array1).astype('float32')
    array2 = numpy.asarray(array2).astype('float32')

    mat1 = cv.fromarray(array1)
    mat2 = cv.fromarray(array2)

    row_margin = int(crop_rate*mat1.rows)
    column_margin = int(crop_rate*mat1.cols)

    rect = ( column_margin, row_margin, mat2.cols-2*column_margin, mat2.rows-2*row_margin )

    mat2_ = cv.GetSubRect( mat2, rect )
   
    template = cv.CreateMat( mat1.rows - mat2.rows + 2*row_margin+1, mat1.cols - mat2.cols + 2*column_margin+1, cv.CV_32FC1 )

    cv.MatchTemplate( mat1, mat2_, template, cv.CV_TM_CCORR_NORMED )

    if ret!='max': return numpy.asarray(template)

    ( minVal, maxVal, minLoc, maxLoc ) = cv.MinMaxLoc(template)

    t = numpy.array(maxLoc) - numpy.array((column_margin,row_margin))

    return numpy.roll(t,1)


def Remap( array_src, array_mapx, array_mapy ):
    '''Do a remapping on numpy arrays'''
    
    flags = cv.CV_INTER_CUBIC + cv.CV_WARP_FILL_OUTLIERS
    fillval = cv.ScalarAll(0)
    
    src = array2cv(array_src.astype('float32'))
    mapx = array2cv(array_mapx.astype('float32'))
    mapy = array2cv(array_mapy.astype('float32'))
    
    shape = array_src.shape
    dest = array2cv(numpy.empty(shape,dtype='float32'))
    
    cv.Remap(src,dest,mapx,mapy,flags,fillval)
    
    array_dest = cv2array(dest)
    
    return numpy.squeeze(array_dest)


def WarpAffine( array_src, transformation2D ):
    '''Do a remapping on numpy arrays'''
    
    flags = cv.CV_INTER_CUBIC + cv.CV_WARP_FILL_OUTLIERS
    fillval = cv.ScalarAll(0)
    
    src = array2cv(array_src.astype('float32'))
    
    shape = array_src.shape
    dest = array2cv(numpy.empty(shape,dtype='float32'))
    
    M = transformation2D.M
    T = transformation2D.T
    
    tf = cv.CreateMat(2, 3, cv.CV_32FC1)
    
    tf[0,0] = M[1,1]
    tf[0,1] = M[1,0]
    tf[0,2] = T[1]
    tf[1,0] = M[0,1]
    tf[1,1] = M[0,0]
    tf[1,2] = T[0]
    
    cv.WarpAffine(src,dest,tf,flags,fillval)
    
    array_dest = cv2array(dest)
    
    return numpy.squeeze(array_dest)


def warp( u, trf2D ):
    '''u is a numpy array to be warped according to a transformation "trf2D"'''

    u = numpy.asarray(u)

    M = trf2D.M
    T = trf2D.T

    nx,ny = u.shape

    warp = cv.CreateMat(2,3,cv.CV_32FC1)

    warp[0,0] = M[1,1]
    warp[0,1] = M[1,0]
    warp[1,0] = M[0,1]
    warp[1,1] = M[0,0]

    warp[0,2] = T[1]
    warp[1,2] = T[0]

    src = cv.fromarray(u.astype('float32'))
    dest = cv.CreateMat(nx,ny,cv.GetElemType(src))

    flags = cv.CV_INTER_CUBIC + cv.CV_WARP_FILL_OUTLIERS

    cv.WarpAffine(src, dest, warp,flags=flags)

    u = numpy.asarray(dest).astype('float')
    u = numpy.where(numpy.isfinite(u),u,0.0)

    return u

    
    
def Resize( array_src, new_shape ):
    '''Resize a numpy array source to dest. Arrays can be conplex.'''
    
    if numpy.iscomplexobj(array_src):
        
        src_r = array2cv(array_src.real.astype('float32'))
        src_i = array2cv(array_src.imag.astype('float32'))
        
        dest_r = array2cv(numpy.empty(new_shape,dtype='float32'))
        dest_i = array2cv(numpy.empty(new_shape,dtype='float32'))
        
        cv.Resize( src_r, dest_r, cv.CV_INTER_CUBIC )
        cv.Resize( src_i, dest_i, cv.CV_INTER_CUBIC )
            
        array_dest = cv2array(dest_r) + cv2array(dest_i)*complex(0.0,1.0)
        
        return numpy.squeeze(array_dest)
    
    else:
    
        src = array2cv(array_src.astype('float32'))
        dest = array2cv(numpy.empty(new_shape,dtype='float32'))
        
        cv.Resize( src, dest, cv.CV_INTER_CUBIC )
            
        array_dest = cv2array(dest)
        
        return numpy.squeeze(array_dest)
