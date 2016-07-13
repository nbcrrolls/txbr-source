"""
.. `modelutil-label`:

The module *modelutil* contains general functions to quickly create/load *3dmod*
models from/to numpy arrays without having to directly use the classes Model, Object
and Contours defined in the module *model*.


"""

import mrc

from model import *

def saveAs3DModel( filename, points, nx=4000, ny=4000, nz=61, xmin=0.0, ymin=0.0, zmin=0.0,
                   sx=1.0, sy=1.0, sz=1.0, doShow=True, mrcFileName=None ):
    '''
    Save a collection of points into a model file called *filename*.

    :param filename: The name of the output 3dmod model file.
    :param points: Can be either a numpy array of shape (n,3), a list of numpy.arrays of similar shape or a list of list of numpy arrays of similar shape
    :param  nx, ny, nz: The size of the model (integer arguments)
    :param xmin, ymin, zmin: The lower boundaryes of the model (float arguments)
    :param sx, sy, sz: The scale of the model (float arguments)
    :param doShow: boolean Open the model into a 3d viewer if True (boolean argument)
    :param mrcFile: If this file is different from None, all the parameters of the model 
        (nx,ny,nz,xmin,ymin,zmin,sx,sy,sz) will match the ones of this MRC file. If doShow
        is selected and mrcFileName is different from None, this file will be opened as well in
        the viewer
    '''

    m = Model( filename )

    if mrcFileName!=None and os.path.exists(mrcFileName):
        f = mrc.MRCFile(mrcFileName)
        m.nx = f.nx
        m.ny = f.ny
        m.nz = f.nz
        m.imgref.xscale_new = f.sx
        m.imgref.yscale_new = f.sy
        m.imgref.zscale_new = f.sz
        m.xoffset = f.x_origin
        m.yoffset = f.y_origin
        m.zoffset = f.z_origin
        m.xmax = f.x_origin + f.xlen
        m.ymax = f.y_origin + f.ylen
        m.zmax = f.z_origin + f.zlen
    else:
        m.nx = nx
        m.ny = ny
        m.nz = nz
        m.imgref.xscale_new = sx
        m.imgref.yscale_new = sy
        m.imgref.zscale_new = sz
        m.xoffset = xmin
        m.yoffset = ymin
        m.zoffset = zmin
        m.xmax = xmin + nx*sx
        m.ymax = ymin + ny*sy
        m.zmax = zmin + nz*sz

    if type(points) is numpy.ndarray:
        obj = m.addNewObject()
        for pt in points:
            cont = obj.addNewContour()
            cont.addPoint( pt[0],pt[1],pt[2] )

    if type(points) is list:
        for item in points:
            obj = m.addNewObject()
            if type(item) is list:
                for _ct in item:
                    cont = obj.addNewContour()
#                    print len(item)
#                    print cont
                    for pt in _ct:
                        cont.addPoint( pt[0],pt[1],pt[2] )
            if type(item) is numpy.ndarray:
                if len(item)!=2:
                    nitems = int(item.size/3)
                    numpy.resize(item,(nitems,3))
                for pt in item:
                    cont = obj.addNewContour()
                    cont.addPoint( pt[0],pt[1],pt[2] )

    m.save()

    if doShow: m.show( withMRCFile=mrcFileName )

    return
    

def loadAllPoints( filename ):
    '''Load all points from a *3dmod* model file in a *numpy* array.
    
    All the points from any contours belonging to any objects of a *3dmod* model file from a 
    file called *filename* are indistinctly grouped in a single *numpy* array of shape (n,3)

    :variable filename: The name of the 3dmod model to scan.
    :returns: A numpy array of shape (n,3) containing all the points from the 3dmod model file.
    '''
    
    m = Model(filename)
    m.loadFromFile()
    
    return m.points()


def loadPoints( filename ):
    '''Load all points from a *3dmod* model file as a list of a list of *numpy* array.

    All the points from any contours belonging to any objects of a *3dmod* model file from a
    file called *filename* are grouped in a list of a list of single *numpy* array of shape (n,3).
    The first list indices the objects, while the sub-list indices contours.

    :variable filename: The name of the 3dmod model to scan.
    :returns: A list of a list of numpy array of shape (n,3) containing the points from
        the 3dmod model file.
    '''

    points = []

    m = Model(filename)
    m.loadFromFile()

    for obj in m.objects:
        contours = []
        for cont in obj.contours:
            contours.append(cont.points)
        points.append(contours)

    return points


def sliceModel( file, u=[0.0,0.0,1.0], X0=[0,0,0], l=[100,100], xlim=None, ylim=None, zlim=None ):
    '''This routine allows to select points of a 3dmod model file that belong
    to a specific slice, without making distinctions of the contour/object they belong
    too. Resulting points are stored in a new model file with an extension ".sliced", each
    point being a single object.

    :param file: The name of the input model file.
    :param u: The normal vector to the slice (a list of 3 float).
    :param X0: A 3D point being a reference for the slice position  (a list of 3 float).
    :param l: A list of 2 float which specify the distance between X0 and the two
        sides of the slices. The sum of the two floats represents the total thickness
        of the slice.
    :param xlim, ylim, zlim: Three 2-tuple that specify boundaries in x, y, z. If None,
        no boundary is taken into account.
    :returns: A numpy array of shape (n,3) containing the positions of the point-beads within the slice.

    '''
    
    template_model = Model(file)
    template_model.loadFromFile()
    XYZ = template_model.points()
    
    output = file + ".sliced"
    
    m = Model(output, template=template_model)
    m.removeAllObjects()
    
    XYZ = sliceVol( XYZ, u=u, X0=X0, l=l, xlim=xlim, ylim=ylim, zlim=zlim )
    
    for pt in XYZ:
        o = m.addNewObject()
        c = o.addNewContour()
        c.addPoints(pt)
        
    m.save()
    
    return XYZ


def sliceVol( XYZ, u=[0.0,0.0,1.0], X0=[0,0,0], l=[100,100], xlim=None, ylim=None, zlim=None):
    '''This routine allows to select points of a collection of 3D points that belong
    to a specific slice, without making distinctions of the contour/object they belong
    too. Resulting points are stored in a new model file with an extension ".sliced", each
    point being a single object.

    :param file: The name of the input model file.
    :param u: The normal vector to the slice (a list of 3 float).
    :param X0: A 3D point being a reference for the slice position  (a list of 3 float).
    :param l: A list of 2 float which specify the distance between X0 and the two
        sides of the slices. The sum of the two floats represents the total thickness
        of the slice.
    :param xlim, ylim, zlim: Three 2-tuple that specify boundaries in x, y, z. If None,
        no boundary is taken into account.
    :returns: A numpy array of shape (n,3) containing the positions of the point-beads within the slice.


    '''

    XYZ = cropModel(XYZ, xlim=xlim, ylim=ylim, zlim=zlim)

    u = numpy.asarray(u)
    X0 = numpy.asarray(X0)
    l = numpy.asarray(l)

    cote = numpy.dot(XYZ-X0,u)

    indices, = numpy.where((cote[:]>=-l[0])*(cote[:]<l[1]))

    return XYZ[indices]


def cropModel( XYZ, xlim=None, ylim=None, zlim=None):
    '''Select points from a numpy array XYZ (of shape (n,3)) within a given box
    
    :param XYZ: A numpy array of shape (n,3) containing the 3D point positions.
    :param xlim, ylim, zlim: Three 2-tuple that specify the boundaries in within the x, y, z
        directions [xlim = (xmin,xmax), ylim = (ymin,ymax) and zlim = (zmin,zmax)].
        If None, no boundary along the axis is taken into account.
    :returns: A copy of a trimmed XYZ without the points outside the box.
    '''

    points = XYZ

    if xlim!=None:
        indices, = numpy.where((points[:,0]>=xlim[0])*(points[:,0]<xlim[1]))
        points = points[indices]

    if ylim!=None:
        indices, = numpy.where((points[:,1]>=ylim[0])*(points[:,1]<ylim[1]))
        points = points[indices]

    if zlim!=None:
        indices, = numpy.where((points[:,2]>=zlim[0])*(points[:,2]<zlim[1]))
        points = points[indices]

    return points