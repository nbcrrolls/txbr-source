import sys,os,os.path,math,logging
import numpy,numpy.linalg,pylab
import scipy,scipy.misc,scipy.interpolate
import mrc,modl,util
import txbr.utilities

WIDTH_REF, HEIGHT_REF = 500, 500

log = logging.getLogger()

class InterfaceSet:
    

    def __init__(self,nx,ny,nz):

        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.constraints = []


    def addConstraint(self,constraint):

        self.constraints.append(constraint)


    def listInterfaces(self):

        lst = numpy.array([constraint.z2() for constraint in self.constraints])
        lst = numpy.sort(lst)
        lst = numpy.unique(lst)

        return lst


    def numberOfInterfaces(self):

        return len(self.listInterfaces())


    def segments(self):

        lst = self.listInterfaces()
        n = len(lst) + 1

        segments = []

        for i in range(n):
            if i==0 and lst[i]!=0: segments.append([0,int(lst[i])-1])
            elif i==0 and lst[i]==0: segments.append([0,0])
            elif i==n-1: segments.append([int(lst[n-2]),int(self.nz-1)])
            else: segments.append([int(lst[i-1]),int(lst[i]-1)])

        return segments


    def listConstraints(self,z2=None):

        if z2!=None:
            lst = [index for index,constraint in enumerate(self.constraints) if constraint.z2()==z2]
        else:
            lst = self.constraints

        return lst


    def empty(self):

        enabled_constraints = [constraint for constraint in self.constraints if constraint.enabled()]

        return len(enabled_constraints)==0


    def numberOfConstraints(self,z2=None):

        return len(self.listConstraints(z2=z2))


    def constraintsMap(self,z2):

        u_ = numpy.row_stack([constraint.u.T for constraint in self.constraints if constraint.z2()==z2])
        v_ = -numpy.row_stack([constraint.v.T for constraint in self.constraints if constraint.z2()==z2])    # Sign convention for the ocv remap

        pos = u_[:,0]*complex(1,0) + u_[:,1]*complex(0,1)
        delta = v_[:,0]*complex(1,0) + v_[:,1]*complex(0,1)

        return (pos,delta)


    def tracks(self,z2):

        in_tracks = numpy.row_stack([constraint.sectionsAt(z2) for constraint in self.constraints if constraint.z2()==z2])
        out_tracks = numpy.row_stack([constraint.sectionsAt(z2-1) for constraint in self.constraints if constraint.z2()==z2])

        return (in_tracks,out_tracks)


class Constraint:

    def __init__(self,contours,order):
        '''Sections should be a list of contour points for every section. Contour points are stored
        in a numpy 2D array, whose second dimension should be 3. Variable order is
        the order we use in the polynomial description of the contours.'''

        self.ncts = len(contours)
        self.sections = [contour.points.copy() for contour in contours]
        self.z = []

        for index in range(self.ncts):
            cotes = numpy.unique(self.sections[index][:,2])
            self.z.append(cotes[0])

        # Make sure the contour are sorted according to their z value

        self.z = numpy.asarray(self.z)

        indices = numpy.argsort(self.z)

        self.z = self.z[indices].tolist()
        self.sections = [self.sections[index] for index in indices]

        self.index1 = int((len(self.sections)-1)/2)
        self.index2 = self.index1 + 1

        self.px,self.py = [],[]

        for contour in contours:
            (px_,py_) = contour.polynomialMapping(order)
            self.px.append(px_)
            self.py.append(py_)

        self.x = []
        self.y = []

        npoints = 100
        t = numpy.linspace(0.0,1.0,num=npoints)

        for index in range(self.ncts):
            self.x.append(self.px[index].eval(t))
            self.y.append(self.py[index].eval(t))

        self.u = (self.x[self.index2],self.y[self.index2])

        if self.index1==0:
            self.v = (self.x[self.index1]-self.x[self.index2],self.y[self.index1]-self.y[self.index2])
        else:
            self.v = (2*self.x[self.index1]-self.x[self.index1-1]-self.x[self.index2],\
                      2*self.y[self.index1]-self.y[self.index1-1]-self.y[self.index2])

        self.u = numpy.asarray(self.u)
        self.v = numpy.asarray(self.v)


    def enabled(self):

        enabled = True

        for index,z_ in enumerate(self.z):
            enabled = enabled and z_ == self.z[0] + index

        return enabled


    def z1(self):

        return self.z[self.index1]


    def z2(self):

        return self.z[self.index2]


    def sectionsAt(self,z):

        index = self.z.index(z)

        return self.sections[index]



def transpose(mapx,mapy):

    map = mapx*complex(1,0) + mapy*complex(0,1)

    return map.transpose().copy()


def decimateSegmentsByTwo(segments):

    new_segments = []

    for segment in segments:
        z0,z1 = segment
        delta = z1-z0
        if delta==1:
            new_segments.appends([z0,z1])
        else:
            mid = int(z0 + delta/2)
            new_segments.appends([z0,mid])
            new_segments.appends([mid+1,z1])

    return new_segments


def dewarp( mrcfileName, modelName, order=4, interface=None, test=False, doPlot=False):
    '''This routine applies a 2D in-plane de-warping.
    mrcfileName: path to the raw stacked file
    modelName: path of the file that contains the constraint contours
    order: polynomial order with which we describe the input contours
    '''

    # Load the reconstruction and the interface model file

    f = mrc.MRCFile(mrcfileName)

    model = modl.Model(modelName)
    model.loadFromFile()

    # Define the map for a no-warp case

    x = numpy.arange(f.nx,dtype='float32')
    y = numpy.arange(f.ny,dtype='float32')

    map0 = transpose(*pylab.meshgrid(x,y))

    # Load informations from the models

    npoints = 100

    interset = InterfaceSet(model.xmax,model.ymax,model.zmax)

    for object in model.objects:

        log.info('Object %i of %s has %i sets of contours!' %(object.indexOfObject,os.path.basename(modelName),object.numberOfContours()))

        constraint = Constraint(object.contours,order)

        if constraint.enabled():
            interset.addConstraint(Constraint(object.contours,order))
        else:
            log.info('Object %i of %s is ignored!' %(object.indexOfObject,os.path.basename(modelName)))


    if interset.empty():

        log.info('No interface for the warping in model %s!' %os.path.basename(modelName))

        return


    log.info('%i different contours to warp on %i interface(s)!' %(interset.numberOfConstraints(),interset.numberOfInterfaces()))
    log.info(interset.listInterfaces())

    for z in interset.listInterfaces():

        log.info('Process Interface: %s' %z)

        basename,extension = os.path.splitext(os.path.basename(mrcfileName))
        map_name = '%s_%i_%i.map' %(basename,z-1,z)

        basename,extension = os.path.splitext(os.path.basename(modelName))
        model_name = '%s_%i_%i%s' %(basename,z-1,z,extension)

        pos,delta = interset.constraintsMap(z)

        init = None
        if os.path.exists(map_name):
            log.info('Load map data from: %s' %map_name)
            mapFile = mrc.MRCFile(map_name)
            init = mapFile.getZSliceAt(0)

        map = strainMap2D( pos, delta, map0, init=init)

        if doPlot:
            in_tracks,out_tracks = interset.tracks(z)
            txbr.utilities.plot_contours(pos,delta,in_tracks,out_tracks,f.nx,f.ny,npoints)
            txbr.utilities.plot_deviation_field(map0,map)

        # Save the map file

        log.info('Save map data from: %s' %map_name)

        mapFile = mrc.MRCFile(map_name)
        mapFile.setHeader(f.nx,f.ny,1,mode=4)
        mapFile.setZSliceAt(0,map)
        mapFile.updateHeader()

    # Do the Remap

    justRemapInterfaces = test    # For testing purpose

    if justRemapInterfaces:

        newfile = remapInterfaces(f,interset)
        model = createModel(modelName,interset)

        log.info('imod %s %s' %(newfile.filename,model.filename))
        os.system('imod %s %s' %(newfile.filename,model.filename))

        log.info('imod %s %s' %(f.filename,model.filename))
        os.system('imod %s %s' %(f.filename,model.filename))

    else:

        newfile = remapFullVolume(f,interset)
        model = createModel(modelName,interset,full=True)

        log.info('imod %s %s' %(newfile.filename,model.filename))
        os.system('imod %s %s' %(newfile.filename,model.filename))


def strainMap2D( pos, delta, map0, init=None):
    '''Calculate the entire displacement field map.
    Variable 'pos' is a complex array that describes the original position points.
    Variables 'delta' is a complex array that describes the position shift.
    '''

    shape_0 = map0.shape
    shape_ref = ( WIDTH_REF, HEIGHT_REF )

    x, y = pos.real*shape_ref[0]/shape_0[0], pos.imag*shape_ref[1]/shape_0[1]
    delta_x, delta_y = delta.real*shape_ref[0]/shape_0[0], delta.imag*shape_ref[1]/shape_0[1]
    
    map0_ = util.Resize(map0, shape_ref)
    
    ux = numpy.column_stack((x,y,delta_x))
    uy = numpy.column_stack((x,y,delta_y))

    if init==None:
        mapx_ = util.mapSurface( shape_ref, ux, init=map0_.real)
        mapy_ = util.mapSurface( shape_ref, uy, init=map0_.imag)
    else:
        init = util.Resize(init, shape_ref)
        u = init - map0_
        mapx_ = util.mapSurface( shape_ref, ux, init=u.real)
        mapy_ = util.mapSurface( shape_ref, uy, init=u.imag)

    map_ = mapx_*complex(1,0) + mapy_*complex(0,1)
    map_ += map0_ 
    
    return util.Resize(map_,shape_0)


def remapFullVolume(f,interfaces):
    '''Remap the MRC file
    '''

    mrcfileName = f.filename

    basename,extension = os.path.splitext(os.path.basename(mrcfileName))
    mrc_name = '%s_remap%s' %(basename,extension)

    log.info('Generate the Full Remap Volume File %s' %mrc_name)

    x = numpy.arange(f.nx,dtype='float32')
    y = numpy.arange(f.ny,dtype='float32')

    map0 = transpose(*pylab.meshgrid(x,y))

    newfile = mrc.MRCFile(mrc_name,template=mrcfileName)

    segments = interfaces.segments()
    log.info('segments: %s' %str(segments))
    segments = numpy.asarray(segments)

    mapf = map0

    for segment in segments:

        log.info('segment: %s' %str(segment))

        basename,extension = os.path.splitext(os.path.basename(mrcfileName))
        map_name = '%s_%i_%i.map' %(basename,segment[0]-1,segment[0])

        if os.path.exists(map_name):
            log.info('Read deviation map from File %s' %map_name)
            mapFile = mrc.MRCFile(map_name)
            mapi = mapFile.getZSliceAt(0)
        else:
            mapi = map0

        indexi = segment[0]
        indexf = segment[1]

        for index in range(indexi,indexf+1):

            log.info('Remap slice #%i' %index)

            im1 = f.getZSliceAt(index).astype('float32')

            if indexi==indexf:
                map = mapi
            else:
                map = ((index-indexi)*mapf + (indexf-index)*mapi)/(indexf-indexi)

            mapx = map.real.astype('float32')
            mapy = map.imag.astype('float32')

            im2 = util.Remap(im1,mapy,mapx)
            
            newfile.setZSliceAt(index,im2)

        indexi = indexf

    newfile.updateHeader()

    return newfile


def remapInterfaces(f,interfaces):
    '''Remap the Interface
    '''

    mrcfileName = f.filename

    basename,extension = os.path.splitext(os.path.basename(mrcfileName))
    mrc_name = '%s_remap%s' %(basename,extension)

    log.info('Generate the Remap Volume File %s' %mrc_name)

    x = numpy.arange(f.nx,dtype='float32')
    y = numpy.arange(f.ny,dtype='float32')

    map0 = transpose(*pylab.meshgrid(x,y))

    newfile = mrc.MRCFile(mrc_name,template=mrcfileName)
    newfile.nz = 2*interfaces.numberOfInterfaces()

    for index,z2 in enumerate(interfaces.listInterfaces()):

        map_name = '%s_%i_%i.map' %(basename,z2-1,z2)

        if os.path.exists(map_name):
            log.info('Read deviation map from File %s' %map_name)
            mapFile = mrc.MRCFile(map_name)
            map = mapFile.getZSliceAt(0)
        else:
            map = map0

        im_i = f.getZSliceAt(z2-1).astype('float32')
        im1 = f.getZSliceAt(z2).astype('float32')

        mapx = map.real.astype('float32')
        mapy = map.imag.astype('float32')

        im_f = util.Remap(im1,mapy,mapx)    # Swap the mapping for x and y

        newfile.setZSliceAt(2*index,im_i)
        newfile.setZSliceAt(2*index+1,im_f)

    newfile.updateHeader()

    return newfile


def createModel(modelPath,interfaceSet,correction_test=True,full=False):
    '''This routine recreates a contour model from the file modelPath'''

    # This routine create a model file that will match a volume composed by the sole interfaces

    basename,extension = os.path.splitext(os.path.basename(modelPath))

    model = modl.Model(modelPath)
    model.loadFromFile()

    new_model_name = '%s_remap%s' %(basename,'.mod')
    newModel = modl.Model(new_model_name,template=model)

    newModel.xmax = model.xmax
    newModel.ymax = model.ymax
    newModel.zmax = model.zmax

    for object in model.objects:
        newObject = modl.Object(new_model_name,object.indexOfObject)
        for index,z2 in enumerate(interfaceSet.listInterfaces()):
            z1 = z2 - 1
            for contour in object.contours:
                if numpy.all(contour.points[:,2]==z1):
                    c_ = contour.copy()
                    if not full: c_.points[:,2] = 2*index
                    newObject.addContour(c_)
                if correction_test and numpy.all(contour.points[:,2]==z1):
                    c_ = contour.copy()
                    if not full: c_.points[:,2] = 2*index + 1
                    else:  c_.points[:,2] += 1
                    newObject.addContour(c_)
                if not correction_test and numpy.all(contour.points[:,2]==z2):
                    c_ = contour.copy()
                    if not full: c_.points[:,2] = 2*index + 1
                    newObject.addContour(c_)
        newModel.addObject(newObject)

    newModel.save()

    return newModel


