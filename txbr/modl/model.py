"""
.. `model-label`:

In the module *model*, classes that redefine *3dmod* contours, objects and
models are re-introduced within the context of a python usage.


"""

import os
import os.path
import re
import logging
import logging.config
import numpy
import scipy.optimize
import pylab
import modfile
import util

log = logging.getLogger()

class _ImageReference:
    '''A class to handle scaling change in IMOD.'''

    def __init__(self):

        self.reset('x')
        self.reset('y')
        self.reset('z')
        
    def reset(self, direction):
        
        if direction=='x':
            self.xscale_old = 1.0
            self.xtranslation_old = 0.0
            self.xrotation_old = 0.0
            self.xscale_new = 1.0
            self.xtranslation_new = 0.0
            self.xrotation_new = 0.0
        if direction=='y':
            self.yscale_old = 1.0
            self.ytranslation_old = 0.0
            self.yrotation_old = 0.0
            self.yscale_new = 1.0
            self.ytranslation_new = 0.0
            self.yrotation_new = 0.0
        if direction=='z':
            self.zscale_old = 1.0
            self.ztranslation_old = 0.0
            self.zrotation_old = 0.0
            self.zscale_new = 1.0
            self.ztranslation_new = 0.0
            self.zrotation_new = 0.0

    def swapAxis(self,axes):
        '''Swap axes. Variable axes is a list of 2 axis index in (0,1,2)'''

        axis1,axis2 = axes

        scale_old = self.xscale_old,self.yscale_old,self.zscale_old
        translation_old = self.xtranslation_old,self.ytranslation_old,self.ztranslation_old
        rotation_old = self.xrotation_old,self.yrotation_old,self.zrotation_old

        scale_new = self.xscale_new,self.yscale_new,self.zscale_new
        translation_new = self.xtranslation_new,self.ytranslation_new,self.ztranslation_new
        rotation_new = self.xrotation_new,self.yrotation_new,self.zrotation_new

        scale_old[axis1],scale_old[axis2] = scale_old[axis2],scale_old[axis1]
        translation_old[axis1],translation_old[axis2] = translation_old[axis2],translation_old[axis1]
        rotation_old[axis1],rotation_old[axis2] = rotation_old[axis2],rotation_old[axis1]

        scale_new[axis1],scale_new[axis2] = scale_new[axis2],scale_new[axis1]
        translation_new[axis1],translation_new[axis2] = translation_new[axis2],translation_new[axis1]
        rotation_new[axis1],rotation_new[axis2] = rotation_new[axis2],rotation_new[axis1]

        self.xscale_old,self.yscale_old,self.zscale_old = scale_old
        self.xtranslation_old,self.ytranslation_old,self.ztranslation_old = translation_old
        self.xrotation_old,self.yrotation_old,self.zrotation_old = rotation_old

        self.xscale_new,self.yscale_new,self.zscale_new = scale_new
        self.xtranslation_new,self.ytranslation_new,self.ztranslation_new = translation_new
        self.xrotation_new,self.yrotation_new,self.zrotation_new = rotation_new
        
        
    def __str__(self):
        
        s = "Image Reference:"
        s += "\nTranslation: (%.1f,%.1f,%.1f)" %(xtranslation_new, ytranslation_new, ztranslation_new)
        s += "\nScale: (%.1f,%.1f,%.1f)" %(xscale_new, yscale_new, zscale_new)
        s += "\nRotation: (%.1f,%.1f,%.1f)" %(xrotation_new, yrotation_new, zrotation_new)

        return s


class Model:
    '''A class to describe a model (with objects, contours and points).'''

    def __init__( self, filename, template=None ):

        self.filename = filename

        self.imgref = _ImageReference()
        
        self.xoffset, self.yoffset, self.zoffset = 0,0,0
        self.xmax, self.ymax, self.zmax = 0,0,0
        self.sx, self.sy, self.sz = 1.0, 1.0, 1.0

        if template!=None:
            self.imgref = template.imgref
            self.xoffset,self.yoffset,self.zoffset = template.xoffset,template.yoffset,template.zoffset
            self.xmax,self.ymax,self.zmax = template.xmax,template.ymax,template.zmax

        self.objects = []


    def loadFromFile( self, filePath=None, keepOnlyPointStructures=False ):
        '''Load the model from an imod file'''

        if filePath==None: filePath = self.filename

        if os.path.exists(filePath): modfile.loadModel( self, filePath )
        else: log.info( 'File {0:s} does not exist'.format(filePath) )
            
        self.removeDisabledObjects()
        
        if keepOnlyPointStructures: self.keepOnlyPointStructures()

        # Do not sort contours here...???
        
    
    def updateImageReference(self, scale=(1.0,1.0,1.0), translation=(0.0,0.0,0.0), rotation=(0.0,0.0,0.0)):
                
        for object in self.objects:
            for contour in object.contours:
                X = self.imgref.xtranslation_new + self.imgref.xscale_new*contour.points[:,0]
                Y = self.imgref.ytranslation_new + self.imgref.yscale_new*contour.points[:,1]
                Z = self.imgref.ztranslation_new + self.imgref.zscale_new*contour.points[:,2]
                contour.points[:,0] = (X - translation[0])/scale[0]
                contour.points[:,1] = (Y - translation[1])/scale[1]
                contour.points[:,2] = (Z - translation[2])/scale[2]
                
        self.imgref.xtranslation_new = translation[0]
        self.imgref.ytranslation_new = translation[1]
        self.imgref.ztranslation_new = translation[2]
        
        self.imgref.xscale_new = scale[0]
        self.imgref.yscale_new = scale[1]
        self.imgref.zscale_new = scale[2]
        

    def sortContours(self):

        for object in self.objects: object.sortContours()


    def numberOfObjects(self):
        '''Return the total number of objects'''

        return len(self.objects)


    def addNewObject(self):
        '''Add an object to the model. Object is indexed by its position in
        the self.objects list.'''

        object = Object(self.filename,len(self.objects))
        self.objects.append(object)

        return object


    def addObject(self,object):
        '''Add an object to the model. Object is indexed by its position in
        the self.objects list.'''

        object.filename = self.filename
        object.indexOfObject = len(self.objects)

        self.objects.append(object)


    def getObject(self,index):
        '''Returns the object at a given index in the self.object list.'''

        return self.objects[index]


    def removeObject(self,index):
        '''Remove an object from the model. Remaining objects are reindexed.'''

        del self.objects[index]

        index = 0
        for object in self.objects:
            object.indexOfObject = index
            for contour in object.contours:
                contour.indexOfObject = index
            index += 1
            
    def removeAllObjects(self):
        
        objects = [object for object in self.objects]
        
        for object in objects: 
            self.removeObject(object.indexOfObject)
        
            
    def removeDisabledObjects(self):
        
        disabled_objects = [object for object in self.objects if object.name=='disabled']

        for object in disabled_objects: self.removeObject(object.indexOfObject)
            
            
    def removeAllPointStructures(self):
        
        point_structures = [object for object in self.objects if object.isAPoint()]
        
        for object in point_structures: self.removeObject(object.indexOfObject)
        

    def keepOnlyPointStructures(self):
        
        non_point_structures = [object for object in self.objects if not object.isAPoint()]
        
        for object in non_point_structures: self.removeObject(object.indexOfObject)
        

    def swapAxis(self,axes):
        '''Swap two axes in the model. Variable axes is a list of 2 axis index in (0,1,2)'''

        axis1,axis2 = axes

        offsets = [self.xoffset,self.yoffset,self.zoffset]
        bounds = [self.xmax,self.ymax,self.zmax]

        offsets[axis1],offsets[axis2] = offsets[axis2],offsets[axis1]
        bounds[axis1],bounds[axis2] = bounds[axis2],bounds[axis1]

        self.xoffset,self.yoffset,self.zoffset = offsets
        self.xmax,self.ymax,self.zmax = bounds

        self.imgref.swapAxis(axes)

        for object in self.objects:
            object.swapAxis(axes)


    def points(self):
        '''Returns an numpy array (of shape (n,3)) of all the points contained in this
        model. Every object is included'''

        pts = [object.points() for object in self.objects if not object.isEmpty()]

        if len(pts)>0:
            return numpy.row_stack(pts)
        else:
            return numpy.zeros((0,3))


    def bounds(self):
        '''Returns the boundary limits of the minimal frame that contains all the
        points of any object in this model.'''

        pts = self.points()

        minimuns = numpy.min(pts,axis=0)
        maximums = numpy.max(pts,axis=0)

        return numpy.row_stack((minimuns,maximums))


    def polynomialMapping( self,order, u_t=[1.0,0.0,0.0], doPlot=False ):
        '''This routine performs a polynomial mapping on every object contained
        in the model.
        '''

        polyMap = []

        for index,object in enumerate(self.objects):
            px,py = object.polynomialMapping(order,u_t=u_t)
            indicesOfViews = object.indicesOfViews()
            polyMap.append((px,py,indicesOfViews))

        if doPlot:

            mapping = []

            for iobject,object in enumerate(self.objects):
                if object.isEmpty(): 
                    continue
                pts = object.points()
                ranges = pts.max(axis=0) - pts.min(axis=0)
                n_t = 5
                n_u = int(ranges[2]+1)
                if object.isAPoint():
                    n_t = 1
                    n_u = 1
                elif object.isALine():
                    n_u = 1
                px = polyMap[iobject][0]
                py = polyMap[iobject][1]
                mapping.append((pts,px,py,n_t,n_u))

            limits = [0.0,self.xmax,0.0,self.ymax]
            label =  'Polynomial Mapping for Surface Traces'

            plot_mapping(mapping,limits=limits,title=label)

        return polyMap


    def tracks(self,n3,u_t=[1.0,0.0,0.0],doPlot=False):
        '''Returns an numpy array (of shape (ntilts,nobject,nn3)) about the contour tracking
        data. Each contour is approximated individually by a polynomial function.'''

        for object in self.objects:
            zvalues = object.zExtensionSet()
            iobject = object.indexOfObject
            for contour in object.contours:
                (pc_x,pc_y) = contour.polynomialMapping(n3,u_t=u_t,doPlot=self.doPlot)
                tracks[itilt,iobject,0,:nn3] = pc_x.coeff[:]
                tracks[itilt,iobject,1,:nn3] = pc_y.coeff[:]

        return tracks


    def save(self):
        '''Save the model in the file self.filename.'''

        head,tail = os.path.split(self.filename)
        if not os.path.lexists(head): os.makedirs(head)

        modfile.saveModel(self)


    def info(self):
        '''Some information about the model.'''

        n = len(self.objects)

        info = '(xmax=%i,ymax=%i,zmax=%i)\n' %(self.xmax,self.ymax,self.zmax)
        info += '(xoffset=%i,yoffset=%i,zoffset=%i)\n' %(self.xoffset,self.yoffset,self.zoffset)
        info += 'Number of Objects: %i\n' %n
        for object in self.objects:
            info += 'Object %i: %s\n%s' %(object.indexOfObject,object.name,object.info())

        return info


    def show( self, withMRCFile=None ):
        '''Open the model with the IMOD viewer'''

        if withMRCFile!=None and os.path.exists(withMRCFile):
            os.system( '3dmod %s %s' %( withMRCFile, self.filename) )
        else:
            os.system( '3dmod -V %s' %self.filename )




class Object:
    '''A class to describe an object in the model'''

    def __init__(self,filename,indexOfObject):

        self.name = None
        self.filename = filename
        self.indexOfObject = indexOfObject
        self.contours = []


    def loadFromFile(self):

        if os.path.exists(self.filename):
            numberOfContours = modfile.numberOfContours(self.filename,self.indexOfObject)
            for i in range(numberOfContours):
                contour = Contour(self.filename,self.indexOfObject,i)
                contour.loadFromFile()
                self.contours.append(contour)


    def numberOfContours(self):
        '''Returns the number of contours that defines this object.'''

        return len(self.contours)


    def addNewContour(self):

        contour = Contour(self.filename,self.indexOfObject,len(self.contours))
        self.contours.append(contour)

        return contour


    def addContour( self, contour ):

        if contour==None:
            contour = Contour( self.filename, self.indexOfObject, len(self.contours) )
        else:
            contour.indexOfObject = self.indexOfObject
            contour.indexOfContour = len(self.contours)

        self.contours.append(contour)

        return contour


    def getContour(self,index):

        return self.contours[index]


    def removeContour(self,index):
        '''Remove a contour from an object. Remaining contours are reindexed.'''

        del self.contours[index]

        for index in range(len(self.contours)):
            contour.indexOfContour = index
            
    def isEmpty(self):
        
        for contour in self.contours:
            if not contour.isEmpty(): return False
        
        return len(self.contours)==0


    def sortContours(self):

        self.contours.sort()


    def zExtensionSet(self):

        pts = [contour.points[:,2] for contour in self.contours]
        zvalues = set(numpy.concatenate(pts,axis=0))

        return zvalues


    def swapAxis(self,axes):
        '''Swap axes. Variable axes is a list of 2 axis index (0,1,2)'''

        for contour in self.contours:
            contour.swapAxis(axes)


    def isAPoint(self):

        isAPoint = True
        for contour in self.contours:
            isAPoint = isAPoint and contour.isAPoint()

        return isAPoint


    def isALine(self):

        return re.match('\s*line\s*', self.name)


    def points(self,sort_points=True):
        '''Returns a numpy array containing all the contour points.
        If keyword sort_points is set to True, then points of each contour are sorted
        along their z,y,x values.'''

        if sort_points:
            pts = [contour.points[numpy.lexsort(contour.points.T).T] for contour in self.contours if not contour.isEmpty()]
        else:
            pts = [contour.points for contour in self.contours if not contour.isEmpty()]

        if len(pts)==0:
            return numpy.zeros((0,3))
        else:
            return numpy.row_stack(pts)


    def bounds(self):
        '''Returns the boundary limits of the minimal frame that contains all the
        points of this object.'''

        pts = self.points()

        minimuns = numpy.min(pts,axis=0)
        maximums = numpy.max(pts,axis=0)

        return numpy.row_stack((minimuns,maximums))


    def minimums(self):

        return numpy.min(contour.points,axis=0)


    def maximums(self):

        return numpy.max(contour.points,axis=0)


    def indicesOfViews( self ):
        '''
        Returns a list containing all the different z value for which there is
        a point that belongs to a contour in the object.
        For an IMOD fiducial file, this z value usually correspond to a tilt index.
        '''

        if len(self.contours)==0:
            return [] 

        indices = [contour.points[:,2] for contour in self.contours]
        indices = numpy.concatenate(indices,axis=0)
        indices = numpy.round(indices)
        indices = set(numpy.array(indices,dtype=int))

        return list(indices)


    def polynomialMapping(self,order,u_t=[1.0,0.0,0.0],doPlot=False):
        """Contours of the object are approximated as a polynomial mapping
        of two parameters (t,u), varying in [0,1]x[0,1]. u parametrizes the
        index of the contour, while varying t moves a point along a contour.
        the 2D vector u_t gives an input of the direction to choose for t
        in the plane of the micrograph.
        """

        order_t = order
        order_u = order

        orders = numpy.array([order_t,order_u])

        if self.isAPoint(): orders[0] = 0

        if len(self.contours)==0:
            return (None,None)
        
        points = numpy.row_stack([contour.points for contour in self.contours])

        # Array of (t,u) parametrization

        one_block = True

        if one_block:

            t = points[:,0]*u_t[0] + points[:,1]*u_t[1]
            u = points[:,2]

            parameters = numpy.column_stack((t,u))

            minimums = numpy.min(parameters,axis=0)
            maximums = numpy.max(parameters,axis=0)

            parameters[:,0] = (parameters[:,0]-minimums[0])/(maximums[0]-minimums[0])
            parameters[:,1] = (parameters[:,1]-minimums[1])/(maximums[1]-minimums[1])

        else:

            parameters = []

            for contour in self.contours:
                pts_ = contour.points.copy()
                minimums = numpy.min(pts_,axis=0)
                maximums = numpy.max(pts_,axis=0)
                pts_[:,0] = (pts_[:,0]-minimums[0])/(maximums[0]-minimums[0])
                pts_[:,1] = (pts_[:,1]-minimums[1])/(maximums[1]-minimums[1])
                pts_[:,2] = contour.indexOfContour/(len(self.contours)-1.0)
                t = pts_[:,0]*u_t[0] + pts_[:,1]*u_t[1]
                u = pts_[:,2]
                parameters.append(numpy.column_stack((t,u)))

            parameters = numpy.row_stack(parameters)

        # Polynomial regression

        nn = util.numberOfTerms(orders,dim=2)

        coeff = numpy.zeros((3*nn))

        def residuals(coeff):
            [coeff_x,coeff_y,coeff_tilt] = numpy.array_split(coeff,3)
            res_x = points[:,0]-util.poleval(coeff_x,orders,parameters)
            res_y = points[:,1]-util.poleval(coeff_y,orders,parameters)
            res_tilt = points[:,2]-util.poleval(coeff_tilt,orders,parameters)
            return numpy.concatenate((res_x,res_y,res_tilt))

        lsq = scipy.optimize.leastsq(residuals,coeff)
        [coeff_x,coeff_y,coeff_tilt] = numpy.array_split(lsq[0],3)

        px = util.Poly(coeff_x,orders)
        py = util.Poly(coeff_y,orders)

        if doPlot:

            n_t,n_u = 500,len(self.contours)
            map = (points,px,py,n_t,n_u)
            mapping = [map]
            label =  'Polynomial Mapping for Surface Traces (Object #%s)' %self.indexOfObject

            plot_mapping(mapping,title=label)

        return (px,py)


    def info(self):
        '''Some information about the object.'''

        info = 'Number of Contours: %i\n' %(self.numberOfContours())

        for i in range(self.numberOfContours()):
            info += '%i -> %s' %(i,self.contours[i].info())

        return info


class Contour:
    '''A class to describe a contour in the model. A contour is a collection of points.'''


    def __init__(self,filename,indexOfObject,indexOfContour):

        self.label = None
        self.filename = filename
        self.indexOfObject = indexOfObject
        self.indexOfContour = indexOfContour
        self.points = numpy.zeros((0,3),dtype=float)


    def loadFromFile(self):

        if os.path.exists(self.filename):
            self.points = numpy.zeros((0,3),dtype=float)
            self.points = modfile.loadPoints(self.filename,self.indexOfObject,self.indexOfContour)


    def copy(self):

        contour = Contour(self.filename,-1,-1)
        contour.points = self.points.copy()

        return contour


    def npoints(self):

        return self.points.shape[0]


    def addPoints(self,points):

        self.points = numpy.row_stack([self.points,points])


    def addPoint(self,x,y,z):

        self.points = numpy.row_stack([self.points,numpy.array([x,y,z])])


    def getPoint( self, index ):

        return self.points[index,0:3]


    def removePoint( self, index ):

        self.points = delete(self.points,index,0)


  

    def getPointAtZ( self, z ):

#        close_points = numpy.zeros((0,3),dtype=float);
#        n = self.npoints()
#        for i in numpy.arange(n):
#            if self.points[i,2]==z:
#                close_points = numpy.row_stack([close_points,self.points[i,0:3]])
#
#        return close_points

        index = numpy.where(self.points[:,2]==z)
        return self.points[index]


    def zExtensionSet(self):

        return set(self.points[:,2])


    def swapAxis(self,axes):
        '''Swap axes. Variable axes is a list of 2 axis index (0,1,2)'''

        axis1,axis2 = axes
        self.points[:,[axis1,axis2]] =  self.points[:,[axis2,axis1]]


    def isAPoint(self):

        return self.points.shape[0]<=1


    def isEmpty(self):

        return self.points.shape[0]==0


    def bounds(self):
        '''Returns the boundary limits of the minimal frame that contains all the
        points of this contour.'''

        n = self.npoints()
        for i in numpy.arange(n):
            min = numpy.amin(self.points,axis=0)
            max = numpy.amax(self.points,axis=0)

        return numpy.row_stack((min,max))


    def polynomialMapping( self, order, u_t=None, doPlot=False):
        """In this routine, contours of an object are approximated as a polynomial function
        of one parameter t, varying between [0,1].
        Variable u_t represents the overall direction when varying t in the 3D space.
        Returns a tuple containing the polynomial approximation for the x and y coordinates
        of points in the contour.
        """

        npoints = self.points.shape[0]

        if npoints==0:  # No point
            return (None,None)

        if npoints==1:  # Single point
            px = util.Poly(self.points[:,0],[0])
            py = util.Poly(self.points[:,1],[0])
            return (px,py)

        # For more than one point...

        orders = numpy.array([npoints<=1 and 0 or min(order,npoints-1)])

        points = self.points[:,:2].copy()
        pts = points[:,:2].copy()

        t = numpy.zeros(pts.shape[0])

        if u_t==None:
            t[1:] = (pts[1:,0]-pts[:-1,0])**2 + (pts[1:,1]-pts[:-1,1])**2
            t = numpy.sqrt(t)
            t = numpy.cumsum(t)
        else:
            t = pts[:,0]*u_t[0] + pts[:,1]*u_t[1]
            t = numpy.sqrt(t**2)

        minimums = numpy.min(t)
        maximums = numpy.max(t)

        if maximums-minimums!=0:
            t = (t-minimums)/(maximums-minimums)
        t.resize(t.size,1)

        # Polynomial regression

        nn = util.numberOfTerms(orders,dim=1)

        coeff = numpy.zeros((2*nn))

        def residuals(coeff):
            [coeff_x,coeff_y] = numpy.array_split(coeff,2)
            res_x = points[:,0]-util.poleval(coeff_x,orders,t)
            res_y = points[:,1]-util.poleval(coeff_y,orders,t)
            return numpy.concatenate((res_x,res_y))

        lsq = scipy.optimize.leastsq(residuals,coeff)
        [coeff_x,coeff_y] = numpy.array_split(lsq[0],2)

        px = util.Poly(coeff_x,orders)
        py = util.Poly(coeff_y,orders)

        if doPlot and False:
            pylab.figure()
            t = numpy.linspace(0,1,100)
            t.resize((t.size,1))
            x = px.eval(t)
            y = py.eval(t)
            pylab.plot(x,y)
            pylab.scatter(points[:,0],points[:,1])
            pylab.title('contour \#%i' %self.indexOfContour)
            #pylab.show()

        return (px,py)


    def __cmp__(self, other):

        if not self.isEmpty() and not other.isEmpty():
            return numpy.min(self.points[:,2])-numpy.min(other.points[:,2])
        elif not self.isEmpty():
            return -1
        elif not other.isEmpty():
            return 1
        else:
            return 0


    def __str__(self):

        if len(self.points)!=0:
            return 'Contour #%i at z=%f' %(self.indexOfContour,numpy.min(self.points[:,2]))
        else:
            return 'Contour #%i' %(self.indexOfContour)



    def info(self):

        return 'Number of Points: %i\n' %self.npoints()



# Utilities

def plot_mapping(mapping,n_t=50,n_u=50,limits=None,title=None,filename=None):
    """Plot contours [data_pts] and their polynomial mapping [px,py].
    """

    x,y = [],[]
    x_pts,y_pts = [],[]

    points = []

    for map in mapping:

        data_pts,px,py,n_t,n_u = map
        
        if px==None or py==None:
            continue

        t = numpy.linspace(0,1,num=n_t)
        u = numpy.linspace(0,1,num=n_u)

        points.append(data_pts)

        for u_ in u:

            pts = numpy.column_stack((t,numpy.ones(n_t)*u_))

            if pts.shape[0]==1:
                x_pts.append(px.eval(pts))
                y_pts.append(py.eval(pts))
            else:
                x.append(px.eval(pts))
                y.append(py.eval(pts))

    if len(x_pts)!=0 and len(y_pts)!=0:
        x_pts,y_pts = numpy.column_stack(x_pts),numpy.column_stack(y_pts)
    if len(x)!=0 and len(y)!=0:
        x,y = numpy.column_stack(x),numpy.column_stack(y)

    if len(points)==0:
        return

    points = numpy.row_stack(points)

    fig = pylab.figure()

    pylab.plot(x_pts,y_pts,'r:.', markersize=0.1)
    pylab.plot(x,y,'r:.', markersize=0.1)
    pylab.scatter(points[:,0],points[:,1])

    pylab.xlabel('x',fontsize=18)
    pylab.ylabel('y',fontsize=18)

    if limits!=None:
        pylab.xlim(limits[0],limits[1])
        pylab.ylim(limits[2],limits[3])

    if title!=None:
        pylab.title(title,fontsize=18)

    if filename!=None:
        pylab.savefig(filename,format='png')


