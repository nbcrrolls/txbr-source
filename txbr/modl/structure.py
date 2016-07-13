import os.path
import scipy.optimize,scipy.linalg,pylab
import numpy,numpy.linalg,numpy.random
import model,util

import logging.config

log = logging.getLogger()

global DEBUG, SMOOTHING_ORDER, USE_MAYAVI

DEBUG = True
SMOOTHING_ORDER = 4
USE_MAYAVI = True

'''
A package for describing the reconstructed tracks/contours
'''
class StructureSet:
    '''A StructureSet is a class allowing to define organize structures, which 
    typically are the patches used during the contour alignment, in a single set.
    Unlike a pure model object, it takes into account all the informations gathered 
    by all the acquisition series'''

    def __init__( self, B, xmax, ymax, directory=None, basename=None ):

        self.directory = directory
        self.basename = basename

        self.nb,self.b = B    # projection map

        self.xmin = 1
        self.xmax = xmax
        self.ymin = 1
        self.ymax = ymax
        self.zmin = -max(xmax,ymax)/2.0
        self.zmax = -self.zmin

        self.nnb = util.numberOfTerms(self.nb,dim=3)
        self.powers_b = util.powerOrder(self.nb,dim=3)

        self.ntilt = self.b.size/2/self.nnb
        self.b = numpy.resize(self.b,(2,self.nnb,self.ntilt))

        self.structures = []


    def clear(self):
        '''Clear the set.'''

        for struture in self.structures:
            structure.set = None

        del self.structures[:]


    def addNewStructure(self,structure):
        '''Add a new structure to the set.'''

        self.structures.append(structure)
        structure.set = self
        structure.indexOfStructure = len(self.structures)-1

        return structure


    def numberOfStructures(self):
        '''Return the number of structures in that set'''

        return len(self.structures)


    def saveModel( self, nx, ny, nz, xoffset=0, yoffset=0, zoffset=0, \
                   scalex=1.0, scaley=1.0, scalez=1.0,
                   filename=None, show=True, direct=False ):
        '''Save the structure set in a model file'''
        
        log.info("StructureSet save model: direct=%s" %direct)
        log.info("(nx,ny,nz): (%i,%i,%i)" %(nx,ny,nz))

        if filename==None:
            filename = os.path.join( self.directory, self.basename + '.mod' )
            
        nx = int(nx)
        ny = int(ny)
        nz = int(nz)
            
#        template = os.path.join( self.directory, '..', self.basename + '.fid' )
#        template = model.Model(template)
#        template.loadFromFile()
        
        offset = numpy.array([ xoffset, yoffset, zoffset ])

        structureModel = model.Model(filename, template=None)
        
#        structureModel.imgref.xscale_new = template.imgref.xscale_new
#        structureModel.imgref.yscale_new = template.imgref.yscale_new
#        structureModel.imgref.zscale_new = template.imgref.zscale_new
        
        structureModel.xmax = nx
        structureModel.ymax = ny
        structureModel.zmax = nz
        
        structureModel.imgref.xscale_new = scalex
        structureModel.imgref.yscale_new = scaley
        structureModel.imgref.zscale_new = scalez

#        structureModel.xoffset = - xoffset
#        structureModel.yoffset = - yoffset
#        structureModel.zoffset = - zoffset
        
        for structure in self.structures:
            
            object = structureModel.addNewObject()
            
            if direct and not structure.isAPoint:
                try:
                    points = structure.points_OC
                except AttributeError: # No initailization for this structure
                    continue
            else:
                points = structure.points
            
            points += offset
                            
            npoints = points.shape[0]
            ncontours = points.shape[1]
                        
            for icontour in range(ncontours):
                contour = object.addNewContour()
                for ipoint in range(npoints):
                    X = points[ipoint,icontour,0]
                    Y = points[ipoint,icontour,1]
                    if True or (X>=0 and X<nx and Y>=0 and Y<ny):
                        contour.addPoints(points[ipoint,icontour,0:3])

        structureModel.save()
        
        if show: structureModel.show()
        
        self.savePointStructures()
        
        
    def savePointStructures( self, filename=None ):
        '''Save the (X,Y,Z) positions of point structures in a file'''
        
        if filename==None:
            filename = os.path.join( self.directory, self.basename + '.markers.txt' )

        f = open(filename,"w")
        
        for structure in self.structures:
            if not structure.isAPoint: 
                continue
            (x,y,z) = structure.points[0,0,0:3]
            f.write('%i  %f  %f  %f\n' %(structure.indexOfStructure,x,y,z))
            
        f.close()


    def plot_structures(self, direct=None):
        '''Plot all the structures present in that set.
        If variable direct is set to True, points that are plotted are the ones
        obtained from the initialization process (dual or real space). 
        If variable direct is set to False, points that are plotted
        represent polynomial approximation surface.
        If variable direct is set to None, both points are plotted.
        '''
        
        box = [ self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax ]
        plot_structures( self.structures, box, direct=direct )


    def info(self):

        info = 'Number of different structures: %i' %self.numberOfStructures()
        
        return info


class Structure:
    '''Class to describe a structure patch in 3D. A structure can either be a point, a line or 
    a 2D surface.'''

    global n_data

    n_data = 10

    def __init__(self, A, C, D, indexOfStructure=0, isAPoint=False, isALine=False, hasTiltViews=True):

        self.indexOfStructure = indexOfStructure

        self.isAPoint = isAPoint
        self.isALine = isALine
        
        self.hasTiltViews = hasTiltViews

        self.update(A,C,D)


    def update(self,A,C,D):
        '''This routine re-evaluate a set of points that belong to the structure
        given its polynomial description.'''

        self.na, self.a = A    # integer, numpy array
        self.nc, self.c = C
        self.nd, self.d = D

        self.nna = util.numberOfTerms(self.na,dim=2)
        self.nnc = util.numberOfTerms(self.nc,dim=1)
        self.nnd = util.numberOfTerms(self.nd,dim=1)

        self.powers_a = util.powerOrder(self.na,dim=2)
        self.powers_c = util.powerOrder(self.nc,dim=1)
        self.powers_d = util.powerOrder(self.nd,dim=1)

        ntilt = self.c.size/2/self.nnc  # Total number of tilts
        self.c.resize(2,self.nnc,ntilt)

        # Generate points of the structure from its polynomial mapping

        orders = numpy.array([self.na,self.na])

        coeff_x = self.a[0:self.nna]
        coeff_y = self.a[self.nna:2*self.nna]
        coeff_z = self.a[2*self.nna:3*self.nna]

        pax = util.Poly(coeff_x,orders)
        pay = util.Poly(coeff_y,orders)
        paz = util.Poly(coeff_z,orders)

        #self.points = self.generate_points(pax,pay,paz,t=(0.25,0.65))
        self.points = self.generate_points(pax,pay,paz)
        

    def generate_points( self, px, py, pz, t=(0.0,1.0), u=(0.0,1.0) ):
        """Given a polynomial mapping parametrized by two variables, this function
        generates points of the surface. orders is a a sequence or numpy array of
        two elements that represent the maximal order used in the polynomial
        parametrization of X, Y and Z.
        """

        if self.isAPoint:
            n_data = 1
            (t0,t1) = (0.0,0.0)
            (u0,u1) = (0.0,0.0)
        elif self.isALine:
            n_data = 10
            (t0,t1) = t
            (u0,u1) = (0.0,0.0)
        else:
            n_data = 20
            (t0,t1) = t
            (u0,u1) = u

        t = numpy.linspace(t0,t1,num=n_data)
        u = numpy.linspace(u0,u1,num=n_data)

        u = numpy.repeat(u,n_data)
        t = numpy.tile(t,n_data)

        parameters = numpy.column_stack((t,u))

        X = px.eval(parameters)
        Y = py.eval(parameters)
        Z = pz.eval(parameters)

        XYZ = numpy.column_stack((X,Y,Z))
        XYZ.resize((n_data,n_data,3))

        return XYZ


    def generate_surface_mapping( self, orders, parameters, XYZ ):
        """Given a set of XYZ points and a set of corresponding (2D) parameters, this function generates
        a polynomial surface mapping of the structure.
        """

        nt,nu,ndim = XYZ.shape[:3]

        XYZ.resize((nt*nu,ndim))

        orders = numpy.asarray(orders)
        nn = util.numberOfTerms(orders,dim=2)

        coeff = numpy.zeros((3*nn))

        def residuals(coeff):
            [coeff_X,coeff_Y,coeff_Z] = numpy.array_split(coeff,3)
            res_X = XYZ[:,0]-util.poleval(coeff_X,orders,parameters)
            res_Y = XYZ[:,1]-util.poleval(coeff_Y,orders,parameters)
            res_Z = XYZ[:,2]-util.poleval(coeff_Z,orders,parameters)
            return numpy.concatenate((res_X,res_Y,res_Z))

        lsq = scipy.optimize.leastsq(residuals,coeff)
        (coeff_X,coeff_Y,coeff_Z) = numpy.array_split(lsq[0],3)

        pX = util.Poly(coeff_X,orders)
        pY = util.Poly(coeff_Y,orders)
        pZ = util.Poly(coeff_Z,orders)

        XYZ.resize((nt,nu,ndim))

        return ( pX, pY, pZ )


    def initializeSurfaceMapping( self, px, py, skipViews, u_t, mode='3D', \
                                  directory=None, basename=None):
        """Given a polynomial mapping for a set of tracks (contours on 2D micrograhs)
        parametrized by two parameters (referred as t and u - the latter is related
        to the tilt angle/number), the function initializeSurfaceMapping generates a 
        polynomial description of the structure itself. 
        Points of the surfaces are then evaluated, either in the direct space when possible
        or in the dual space - occluded contours method (Kang 2001). It could be done either 
        in 2D (mode='2D') or 3D (mode='3D'). In 2D, contour points of the surface are 
        generated for different values of t, while in 3D point of the whole surface are 
        generated with the dual space method outlined by Kang. 
        Points of the surface are then stored in the variable self.points_OC
        A least square regression method is then used to obtain the best polynomial 
        mapping (of order [self.na,self.na]) of the surface in 3D. 
        Results are stored in the variables self.a and self.points
        Variable directory: specifies a directory where to store some running time informations
        Variable basename: specifies the basename series for which this structure is built.
        """

        log.info('Init Surface Mapping for structure #%i (u_t=%s)' %(self.indexOfStructure,str(u_t)))

        if mode!='2D': mode='3D'
        
        if numpy.max(px.order)<=0 and numpy.max(py.order)<=0:
            log.info('Cannot initialize structure with polynomials or order 0')
            return
        
        indices = numpy.where(skipViews==False)[0]
        ntilt = indices.size
        
        if ntilt==0:
            log.warning('Surface #%i has no projection view!' %self.indexOfStructure)
        
        filename = os.path.join(directory,basename)

        # 1) First evaluate tangents for the projected contours

        t0,t1 = 0.0,1.0
        u0,u1 = 0.0,1.0
        
        t,u = numpy.meshgrid(numpy.linspace(t0,t1,n_data),numpy.linspace(u0,u1,ntilt))
        parameters = numpy.column_stack((t.ravel(),u.ravel()))

        x0 = px.eval(parameters)
        y0 = py.eval(parameters)

        xy_0 = numpy.column_stack((x0,y0))

        dpxdt = px.der((1,0))
        dpydt = py.der((1,0))

        xprime = dpxdt.eval(parameters)
        yprime = dpydt.eval(parameters)

        tangents_p = numpy.zeros((ntilt*n_data,3))    # Coefficients of the projected tangents

        tangents_p[:,0] = yprime[:]*x0[:] - y0[:]*xprime[:]
        tangents_p[:,1] = -yprime[:]
        tangents_p[:,2] = xprime[:]
        
        tangents_p.resize((ntilt,n_data,3))

        # 2) extension to the volume: evaluation of the tangent planes to the surface patches

        tangents = numpy.zeros((ntilt,n_data,4))

        P = self.set.b[:,0:4,:]

        tangents[:,:,0] += tangents_p[:,:,0]
        for itilt in range(ntilt):
            tangents[itilt,:,:] += numpy.tensordot(tangents_p[itilt,:,1:],P[:,:,indices[itilt]],axes=([1],[0]))        
        
        # 3) Pull the tangents and generate points from the surface

        if mode=='2D':
            
            # We consider the tilt axis is the same for each series here...eventually to correct...
            u_t = numpy.average(numpy.array(u_t),axis=0)
            
            n3 = numpy.array([0.0,0.0,1.0])
            n1 = u_t/numpy.dot(u_t,u_t)
            n2 = numpy.cross(n3,n1)
            
            log.debug('u1=%s  u2=%s  u3=%s' %(str(n1),str(n2),str(n3)))

            P = numpy.row_stack((n1,n2,n3))
            Pinv = numpy.linalg.inv(P)
            
            tangents[:,:,1:] = numpy.tensordot(tangents[:,:,1:],P,axes=([2],[1]))

            XY_0 = numpy.tensordot( xy_0, P[:2,:2], axes=([1],[1]) )
            X0 = numpy.resize(XY_0[:,0],(ntilt,n_data))
            
            tangents[:,:,0] += tangents[:,:,1]*X0
            
            tangents = numpy.delete(tangents,1,2)
    
            # order of the algebraic surface in the space of the tangent coefficients
            if self.isAPoint:
                order = 2
            else:
                order = 4
                
            points2D = self.buildSurface2D( tangents, order, tilts=indices, filename=filename )
            
            s = points2D.shape # shape can be (ntilt-1,n_data,2) (direct evaluation) or (ntilt,n_data,2) (dual space evaluation)
            
            X0 = X0[:s[0]].ravel()
            
            points2D.resize((s[0]*s[1],2))
            
            points = numpy.column_stack((X0,points2D))
            points = numpy.tensordot(points,Pinv,axes=([1],[1]))
            
            points.resize((s[0],s[1],3))
            
        elif mode=='3D':

            order = 4    # order of the algebric surface in the space of the tangent coefficients
            
            points = self.buildSurface3D( tangents, order, tilts=indices, filename=filename)
            
        self.indicesOfTilts = indices
        self.points_OC = points

        # Generate a polynomial mapping of the surface of the surface and points of the mapping

        orders = [self.na,self.na]

        ntilt_eff = points.shape[0]
        
        t,u = numpy.meshgrid(numpy.linspace(t0,t1,n_data),numpy.linspace(u0,u1,ntilt_eff))
        parameters = numpy.column_stack((t.ravel(),u.ravel()))

        [pax,pay,paz] = self.generate_surface_mapping( orders, parameters, points )

        # Store the results in the structure

        self.a = numpy.concatenate((pax.coeff,pay.coeff,paz.coeff))
        self.points = self.generate_points(pax,pay,paz,t=(t0,t1))
        
        
    def buildSurface2D( self, tangents, order, tilts=None, mode='direct', filename=None):
        """2D Dual Space Approach. Occluded contour reconstruction. 
        Variable tangents: is an numpy array containing the tangent line
        coefficients to the curve (2D) to reconstruct. For the volume
        parametrization, the contours correspond to t=cte.
        Variable order: order of the polynomial approximation we use to 
        describe the surface.
        Variable mode: specifies is the estimated structure (curve here in 2D) will 
        be done either in the dual space (space of the tangent coefficients) or in the 
        direct space.
        """
               
        global DEBUG, SMOOTHING_ORDER
               
        ntilt, ndata = tangents.shape[:2]
        smoothing_order = SMOOTHING_ORDER
        
        if tilts==None:
            tilts = numpy.arange(ntilt)
            
        tilts = numpy.resize(tilts,(ntilt,1))
        
        # Smooth out the tangent coefficients over the tilt series with a polynomil regression
        smooth_tangents = smooth(tilts,tangents,smoothing_order)
        
        plot_tangents_coefficients( smooth_tangents, ref=tangents, tilts=tilts, \
                                    patch=self.indexOfStructure, filename=filename )
        
        tangents = smooth_tangents
    
        def plot_contour2D( idata, pts, redraw=False ):
            text=r'Structure \#%i   Section Data \#%i' %(self.indexOfStructure+1,idata+1)
            limits = [ self.set.xmin, self.set.xmax, 2*self.set.zmin, 2*self.set.zmax ]
            plotContour2D( tangents[:,idata,:],pts[:,idata,:], figtext=text, limits=limits, \
                           section=idata+1, patch=self.indexOfStructure+1, redraw=redraw, filename=filename)

        log.info('Build the surface for each section perpendicular to the tilt axis! (mode: %s)' %mode)

        if mode=='direct':

            # In the direct method, a point of the contour is approximated by the
            # intersection of two consecutives tangents.

            pts = numpy.zeros((ntilt-1,n_data,2))

            for idata in range(n_data):
                log.debug('Data section #%i' %idata)
                for itilt in range(ntilt-1):
                    A = tangents[itilt:itilt+2,idata,1:]
                    B = -tangents[itilt:itilt+2,idata,0]
                    #print 'tilt #%i: %s' %(itilt,str(A))
                    try:
                        Ainv = numpy.linalg.inv(A)
                    except numpy.linalg.linalg.LinAlgError:
                        print 'Singular Matrix at tilt: %i' %itilt
                        print A
                    pt = numpy.dot(Ainv,B)
                    pts[itilt,idata,:] = pt[:]

            if DEBUG:
                for idata in range(n_data):
                    plot_contour2D( idata, pts, redraw=idata!=0)

            return pts


        if mode=='dual':
            
            # In the dual method, we used the homogeneous nature of the tangents
            # to delimitate the surfaces

            ndim = 3
            nn_start = util.numberOfTerms(order-1,dim=ndim)
            nn_end = util.numberOfTerms(order,dim=ndim)

            pts = numpy.zeros((ntilt*n_data,2))

            coeff = numpy.zeros((nn_end-nn_start))

            for idata in range(n_data):

                p = util.Poly(numpy.zeros((nn_end)),numpy.array([order]*ndim),min_order=order)

                X = p.X(tangents[:,idata,:])
                V = numpy.tensordot(X,X,axes=([1],[1]))

                evalues,evectors = scipy.linalg.eigh(V)

                log.info('Occluded Contour (2D) #%i: fmin=%e' %(idata,min(evalues)))

                esortvalues = numpy.argsort(evalues)
                coeff = evectors[:,esortvalues[1]]

                #coeff = evectors[:,numpy.argmin(evalues)]

                p.coeff[nn_start:nn_end] = coeff[:]

                p1 = p.der([1,0,0])
                p2 = p.der([0,1,0])
                p3 = p.der([0,0,1])

                der1 = p1.eval(tangents[:,idata,:])
                der2 = p2.eval(tangents[:,idata,:])
                der3 = p3.eval(tangents[:,idata,:])

                pts[idata::n_data,0] = der2[:]/der1[:]
                pts[idata::n_data,1] = der3[:]/der1[:]

            pts.resize((ntilt,n_data,2))

            if DEBUG:
                for idata in range(n_data):
                    plot_contour2D(idata,pts)

            return pts

        return None


    def buildSurface3D( self, tangents, order, tilts=None, filename=None ):
        """Occluded surface reconstruction. tangents is an ndarray of the
        tangent planes to the surface to reconstruct.
        """
        
        global DEBUG, SMOOTHING_ORDER
               
        ntilt, ndata, ndim = tangents.shape
        smoothing_order = SMOOTHING_ORDER
        
        tilts = numpy.resize(numpy.arange(ntilt),(ntilt,1))
        
        # Smooth out the tangent coefficients over the tilt series with a polynomil regression
        smooth_tangents = smooth(tilts,tangents,smoothing_order)
        
        plot_tangents_coefficients( smooth_tangents, ref=tangents, tilts=tilts, \
                                    patch=self.indexOfStructure, filename=filename )
        
        tangents = smooth_tangents
        
        tangents.resize(ntilt*n_data,ndim)

        nn_start = util.numberOfTerms(order-1,dim=ndim)
        nn_end = util.numberOfTerms(order,dim=ndim)

        p = util.Poly(numpy.zeros((nn_end)),numpy.array([order]*ndim),min_order=order)

        X = p.X(tangents)

        V = numpy.tensordot(X,X,axes=([1],[1]))

        evalues,evectors = scipy.linalg.eigh(V)

        coeff = evectors[:,numpy.argmin(evalues)]

        p.coeff[nn_start:nn_end] = coeff[:]

        p1 = p.der([1,0,0,0])
        p2 = p.der([0,1,0,0])
        p3 = p.der([0,0,1,0])
        p4 = p.der([0,0,0,1])

        der1 = p1.eval(tangents)
        der2 = p2.eval(tangents)
        der3 = p3.eval(tangents)
        der4 = p4.eval(tangents)

        XYZ = numpy.zeros((ntilt*n_data,3))

        XYZ[:,0] = der2[:]/der1[:]
        XYZ[:,1] = der3[:]/der1[:]
        XYZ[:,2] = der4[:]/der1[:]

        points = XYZ.reshape((ntilt,n_data,3))

        return points

#--------------------------------Helper Functions ----------------------------------------

def smooth( x, y, order ):
    
    from util import polynomial_regression
    
    n = x.size
    x_ = numpy.resize(x,(n,1))

    y_shape= y.shape
    
    if y_shape[0]!=n: raise ValueError
    
    if len(y_shape)==1: 
        y_ = numpy.resize(y,(n,1))
    else:
        y_ = numpy.resize(y,(n,y.size/n))
    
    order_ = numpy.array([order])
    
    smooth_y = numpy.empty_like(y_)
    
    for i in range(y.size/n):
        p = polynomial_regression(x_,y_[:,i],order_)
        smooth_y[:,i] = p.eval(x_)
    
    smooth_y.resize(y_shape)
    
    return smooth_y

# Plotting Functions

def plot_structures( structures, box, direct=None ):
    '''Plot structures in 3D with MayaVi. Structures are represented by their
    polynomial interpolation.'''
    
    global USE_MAYAVI
    
    if not USE_MAYAVI: return

    if direct:
        log.info('Plot structures in 3D with MayaVi from initial estimates!')
    else:
        log.info('Plot structures in 3D with MayaVi')

    try:
        
        # We separate the point-like and patch structures for performance issue
        
        point_structures = [ structure for structure in structures if structure.isAPoint ]
        patch_structures = [ structure for structure in structures if not structure.isAPoint ]
        
        xmin,ymin,zmin,xmax,ymax,zmax = box
        
        def plotXYZ( XYZ, color ):
            
            index = numpy.where((XYZ[:,:,0]>=xmin) & (XYZ[:,:,0]<=xmax) & \
                                (XYZ[:,:,1]>=ymin) & (XYZ[:,:,1]<=ymax) & \
                                (XYZ[:,:,2]>=zmin) & (XYZ[:,:,2]<=zmax))
            XYZ = XYZ[index]
            
            s = mlab.points3d( numpy.ravel(XYZ[:,0]), numpy.ravel(XYZ[:,1]), numpy.ravel(XYZ[:,2]), \
                                       color=color, scale_factor=10.0 )

        # do the plotting

        from enthought.mayavi import mlab

        fig = mlab.figure()
        
        if direct==None or direct==True:
            
            # Plot the point structures
            
            try:
                XYZ0 = [ structure.points_OC for structure in point_structures ]
            except AttributeError:
                XYZ0 = [ structure.points for structure in point_structures ]
                    
            if len(XYZ0)!=0: plotXYZ( numpy.row_stack(XYZ0), (1.0,0.0,0.0))
            
            # Plot the patch structures
            
            for structure in patch_structures:
                
                try:
                    XYZ0 = structure.points_OC
                except AttributeError:
                    XYZ0 = structure.points
                    
                plotXYZ( XYZ0, (1.0,1.0,0.0)) 
                    
        if direct==None or direct==False:
            
            XYZ = [ structure.points for structure in point_structures ]
            
            if len(XYZ)!=0: plotXYZ( numpy.row_stack(XYZ), (1.0,0.0,0.0))
            
            for structure in patch_structures:
                XYZ = structure.points
                plotXYZ( structure.points, (1.0,0.0,1.0))
                
        mlab.show()

    except ImportError:

        print "Cannot complete 3d plotting:", ImportError
            
            
def plot_tangents_coefficients( tgts, ref=None, tilts=None, patch=None, filename=None ):
    '''Plot the coefficients of tangents that are used to reconstruct a structure from
    its contours.'''
    
    n, n_data,dim = tgts.shape
    
    tilts = numpy.asarray(tilts).ravel()
    
    if tilts==None:
        tilts = numpy.arange(n)
    
    if n_data<=5:
        step = 1
    else :
        step = int(n_data/5)
    
    pylab.figure()

    if dim==3:
        for i_data in range(0,n_data,step):
            pylab.plot(tilts,tgts[:,i_data,0],'r')
            pylab.plot(tilts,tgts[:,i_data,1],'g')
            pylab.plot(tilts,tgts[:,i_data,2],'b')
            if ref!=None:
                pylab.scatter(tilts,ref[:,i_data,0],c=(1.0,0.0,0.0))
                pylab.scatter(tilts,ref[:,i_data,1],c=(0.0,1.0,0.0))
                pylab.scatter(tilts,ref[:,i_data,2],c=(0.0,0.0,1.0))
                
    if dim==4:
        for i_data in range(0,n_data,step):
            pylab.plot(tilts,tgts[:,i_data,0],'r')
            pylab.plot(tilts,tgts[:,i_data,1],'g')
            pylab.plot(tilts,tgts[:,i_data,2],'b')
            pylab.plot(tilts,tgts[:,i_data,3],'m')
            if ref!=None:
                pylab.scatter(tilts,ref[:,i_data,0],c=(1.0,0.0,0.0))
                pylab.scatter(tilts,ref[:,i_data,1],c=(0.0,1.0,0.0))
                pylab.scatter(tilts,ref[:,i_data,2],c=(0.0,0.0,1.0))
                pylab.scatter(tilts,ref[:,i_data,3],c=(1.0,0.0,1.0))
                
    pylab.xlim(numpy.min(tilts),numpy.max(tilts))
    pylab.xlabel("Tilt Index")
    
    if patch==None:
        pylab.title(r"Tangent Coefficients vs Tilts")
    else:
        pylab.title(r"Tangent Coefficients vs Tilts (Structure \#%i)" %patch)
        
    if filename!=None and patch!=None:
        pylab.savefig( '%s_tgts_%i.png' %(filename,patch), format='png' )
        

def plot_2D_tangents(tangents,tangents_2,xmin,xmax):

    pylab.figure()

    ntilt = tangents.shape[0]
    n_data = tangents.shape[1]

    xpts = numpy.linspace(xmin,xmax,num=2)

    for itilt in range(ntilt):
        t0 = tangents[itilt,0]
        tx = tangents[itilt,1]
        ty = tangents[itilt,2]
        T2d = numpy.array([-tx/ty,-t0/ty])
        ypts = numpy.polyval(T2d,xpts)
        pylab.plot(xpts,ypts,'r:')

    ntilt = tangents_2.shape[0]
    n_data = tangents_2.shape[1]

    xpts = numpy.linspace(xmin,xmax,num=2)

    for itilt in range(ntilt):
        t0 = tangents_2[itilt,0]
        tx = tangents_2[itilt,1]
        ty = tangents_2[itilt,2]
        T2d = numpy.array([-tx/ty,-t0/ty])
        ypts = numpy.polyval(T2d,xpts)
        pylab.plot(xpts,ypts,'g:')
        
        
def plotContour2D( tangents, XY=None, patch=None, swapaxes=True, labels=['X','Z'], \
                   limits=[0,1000,-400,400], section=None, figtext=None, \
                   redraw=False, filename=None ):
    """Plot tangents and points for a 2D contour.
    """

    if swapaxes:
        (t0,tx,ty) = numpy.array_split(tangents,3,axis=1)
        tangents = numpy.column_stack((t0,ty,tx))
        limits = [ limits[2], limits[3], limits[0], limits[1] ]
        labels = [ labels[1], labels[0] ]
        
    if XY!=None and swapaxes:
        (x,y) = numpy.array_split(XY,2,axis=1)
        XY = numpy.column_stack((y,x))
    
    if not redraw:
        pylab.figure()
    else:
        pylab.gcf()
        
    ntilt = tangents.shape[0]
    n_data = tangents.shape[1]

    xpts = numpy.linspace(limits[0],limits[1],num=2)

    for itilt in range(ntilt):
        t0 = tangents[itilt,0]
        tx = tangents[itilt,1]
        ty = tangents[itilt,2]
        T2d = numpy.array([-tx/ty,-t0/ty])
        ypts = numpy.polyval(T2d,xpts)
        pylab.plot(xpts,ypts,'r:')

    if XY!=None:
        pylab.plot(XY[:,0],XY[:,1], 'g-o')
        
    center = numpy.mean( XY, axis = 0 )
    
    half_width = (limits[1] - limits[0])/2.0
    half_height = (limits[3] - limits[2])/2.0
    
    limits = [ center[0] - half_width, \
               center[0] + half_width, \
               center[1] - half_height, \
               center[1] + half_height ]

    pylab.axis(limits)

    pylab.xlabel( labels[0], fontsize=28 )
    pylab.ylabel( labels[1], fontsize=28 )

    if figtext!=None:
        pylab.title(figtext,fontsize=28)

    if filename!=None and patch!=None and section!=None:
        pylab.savefig( '%s_tgts_%i-%i.png' %(filename,patch,section), format='png' )

