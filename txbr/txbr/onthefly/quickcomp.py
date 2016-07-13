"""
.. `quickcomp-label`:

The module *quickcomp* allows to map the corresponding beads between tilt series.


"""

import os
import sys
import collections
import numpy
import scipy.spatial
import scipy.optimize
import pylab
import numpy.fft
import numpy.linalg
import scipy.ndimage.interpolation
import scipy.ndimage.filters
import mrc
import modl
import util
import txbr.onthefly.phasecorr
import txbr.onthefly.otfutil

from txbr import log
from otfutil import changeToPolarCoordinate

CACHE_FILE = "reltrans.txt"

def __loadTransformations__( project, raw=False ):
    """Load the relative transformations between series within a project"""

    log.info("Load Transformations")

    cache = os.path.join( project.multser_dir, CACHE_FILE )

    if os.path.lexists(cache):
        file = open(cache)
        transformations = eval(file.read())
    else:
        transformations = {}

    if raw: return transformations

    # Otherwise, modify the transformations

    transformations_ = {}

    for label,tranformation in transformations.iteritems():
        ext_a,ext_b = label.split("-")
        indexOfSeries_a = project.indexOfSeries(ext_a,fromext=True)
        indexOfSeries_b = project.indexOfSeries(ext_b,fromext=True)
        log.info("Extensions: [{0:s},{1:s}] -> Indices: [{2:},{3:}]".format(ext_a,ext_b,indexOfSeries_a,indexOfSeries_b))
        if (indexOfSeries_a==-1 or indexOfSeries_b==-1): continue
        label = "{0:}-{1:}".format(indexOfSeries_a,indexOfSeries_b)
        transformations_[label] = util.AffineTransformation3D(tranformation)

    return transformations_


def __saveTransformations__( project, transformations, load=True ):
    """Save the relative transormations between series of a project"""

    log.info("Save Transformations")

    transformations_ = load and __loadTransformations__(project,raw=True) or {}

    cache = os.path.join( project.multser_dir, CACHE_FILE )

    for label,tranformation in transformations.iteritems():
        indexOfSeries_a,indexOfSeries_b = label.split("-")
        ext_a = project.extensions[int(indexOfSeries_a)]
        ext_b = project.extensions[int(indexOfSeries_b)]
        label_ = "{0:s}-{1:s}".format(ext_a,ext_b)
        M =  numpy.column_stack((tranformation.T,tranformation.M)).ravel().tolist()
        transformations_[label_] = M

    file = open(cache, 'w')
    file.write(str(transformations_))
    file.close()


def scanrot( uref, u, angles=(0,2.0*numpy.pi,360), zoom=1.0, verbose=True,
             label=None, work_directory=None, display_slices=False, display_peak=False ):
    """Find the optimal corresponding rotation between two images [uref,u]. The procedure
    is done by cross-correlating the reference image *uref* and a set of rotated images
    for *u*. The best rotation should correspond to a cross-correlating peak.

    :param uref: A 2D numpy array. The reference image.
    :param u: A 2D numpy array. The image to compare with.
    :param angles: A list of rotation angles between the two images to investigate.
    :param zoom: The procedure can be performed on reduced images if *zoom*<1.

    """

    uref = (uref-numpy.mean(uref))/numpy.std(uref)
    u = (u-numpy.mean(u))/numpy.std(u)

    center = numpy.array(uref.shape)/2.0

    if numpy.any(zoom!=1.0):
        uref = scipy.ndimage.interpolation.zoom( uref, zoom )
        u = scipy.ndimage.interpolation.zoom( u, zoom )

    if uref.shape!=u.shape:
        nx = min(uref.shape[0],u.shape[0])
        ny = min(uref.shape[1],u.shape[1])
    else:
        nx,ny = uref.shape

    theta = numpy.linspace(*angles)
    
    ox = int((numpy.sqrt(2.0)-1.0)/2.0/numpy.sqrt(2.0)*nx)
    oy = int((numpy.sqrt(2.0)-1.0)/2.0/numpy.sqrt(2.0)*ny)
    
    lx = nx-2*ox
    ly = ny-2*oy

    stds = []

    peak0 = 0.0

    for index,angle in enumerate(theta):
        
        v = scipy.ndimage.interpolation.rotate(u,numpy.degrees(angle))
        
        o1,o2 = numpy.array(v.shape)/2.0
        x1 = o1-lx/2
        x2 = x1 + lx
        y1 = o2-ly/2
        y2 = y1 + ly
        
        u1 = uref[ox:-ox,oy:-oy]
        u2 = v[x1:x2,y1:y2]

        phase = txbr.onthefly.phasecorr.getPhaseCorrelation( u1, u2, apply_med_filt=None )
        ( tx,ty,std,peak ) = txbr.onthefly.phasecorr.maximum( phase, mode='regular')

        if verbose:
            log.info("#{0:3}   Angle: ({1:5.1f}) {2:5.1f} -> {3:+6.2f}  {4:+6.2f}    {5:.2e}  {6:.3f}".format(index,angle,numpy.degrees(angle),tx,ty,std,peak))
           
        if peak>=peak0:
            peak0 = peak
            tx0 = tx
            ty0 = ty
            angle0 = angle

        stds.append((numpy.degrees(angle),peak))

    log.info("Final Angle: {0:.2f}   shift [{1:.1f},{2:.1f}]".format(numpy.degrees(angle0),tx0,ty0))

    # Final Rotation parameter

    M = numpy.array([[numpy.cos(angle0),-numpy.sin(angle0)], [numpy.sin(angle0),numpy.cos(angle0)]])
    T =  center  - numpy.dot(M,center) + [ tx0/zoom , ty0/zoom ]

    rot = util.AffineTransformation2D( numpy.column_stack((T,M)))

    if not os.path.lexists(work_directory): os.makedirs(work_directory)

    # Reference slices comparison

    c = numpy.array(u.shape)/2.0
    offset = c - numpy.dot(M.T,c) - numpy.dot(M.T,(tx0,ty0))
    v = scipy.ndimage.interpolation.affine_transform( u, M.T, offset=offset)

    filename = os.path.join( work_directory, "slices_{0:s}.mrc".format(label) )

    nx,ny = uref.shape

    f = mrc.MRCFile(filename)
    f.setHeader(nx,ny,2)
    f.setZSliceAt(0,v)
    f.setZSliceAt(1,uref)
    f.updateHeader()

    if display_slices: os.system("imod {0:s}".format(filename))

    # Plot the correlation peak

    peak_png = os.path.join( work_directory, "peak_{0:s}.png".format(label) )
    stds = numpy.asarray(stds)
    pylab.figure()
    pylab.plot(stds[:,0],stds[:,1])
    pylab.xlabel(r"$\theta$")
    pylab.ylabel(r"Correlation Peak")
    pylab.title(r"Rotation angle {0:s}: {1:.1f}$^\circ$".format(label,numpy.degrees(angle0)))
    pylab.savefig(peak_png)
    
    if display_peak: pylab.show()

    return rot


def evaluateAngle( u1, origin1, u2, origin2 ):
    '''Evaluate directly the rotation angle between two images (2d ndarray) by
    cross-correlating the images in polar coordinates. Two corresponding points
    *origin1* and *origin2* need to be known rather accurately for the method
    to be efficient. If not, use the routine *scanrot*.

    :param u1: First image as a 2d numpy array.
    :param origin1: Reference point in the first image.
    :param u2: First image as a 2d numpy array.
    :param origin2: Reference point in the second image.

    '''

    origin1 = numpy.squeeze(origin1)
    origin2 = numpy.squeeze(origin2)

    v1 = changeToPolarCoordinate(u1,origin=origin1)
    v2 = changeToPolarCoordinate(u2,origin=origin2)

    tx, ty, widthx, widthy, theta, a, b =  txbr.onthefly.phasecorr.correlate( v1, v2, crop_rate=0 )

    return tx*2.0*numpy.pi/1000


def minimizeE( pts1, pts2, M0 ):
    """
    Attempt to find the correspondence between two set of points by minimizing
    a free energy like quantity.

    :param pts1: The first set of points.
    :param pts1: The second set of points.

    """


    def f( x ):

        M = util.AffineTransformation2D(numpy.resize(x.copy(),(2,3)))

        pts1_ = pts1
        pts2_ = M.forward(pts2)

        d = scipy.spatial.distance.cdist(pts2_,pts1_)

        e = 1.0 - numpy.exp(-(d/5.0)**2)
        e = numpy.mean(e)

        log.info("tx,ty={0:.2f},{1:.2f}  M=[[{2:.2f},{3:.2f}],[{4:.2f},{5:.2f}]]   {6:.3f}".format(M.T[0],M.T[1],M.M[0,0],M.M[0,1],M.M[1,0],M.M[1,1],e))

        return e


    x0 = [ M0[0,0], M0[0,1], M0[0,2], M0[1,0], M0[1,1], M0[1,2] ]
    
    xopt = scipy.optimize.fmin_cg( f,x0, gtol=1e-7, epsilon=1e-10 )

    return numpy.resize(xopt,(2,3))


def rotateModel( model, rot2d, output ):
    """Rotate all the points of a model in the plane of their two first coordinates.
    The output model is saved in the file specified by the parameter *output*. If nothing is
    specified the character '_' is appended to the original model name
    """

    m = modl.Model(model)
    m.loadFromFile()

    for obj in m.objects:
        for cont in obj.contours:
            cont.points[:,:2] = rot2d.forward(cont.points[:,:2].copy())

    if output!=None:
        m.filename = output
    else:
        m.filename = m.filename + "_"

    m.save()


def compare3DModel( model1, model2, rot3D, dm=2.5, miniter=10, maxiter=50, label=None,
                    work_directory=None, displayModel=False, displayResidual=False ):
    """Compare the 3D locations of two sets of markers. Markers are stored in two
    differnt models *model1* and *model2*, and the 3D transformation *rot3D* is an
    estimation of the relative transformation that map the first set to the second set.
    This mapping and the corresponding transformation are refined through an iterative
    process (with at most *maxiter* iterations). In this process, two points are considered
    identical if the final distance is less than *dm*.
    """

    m1 = modl.Model(model1)
    m1.loadFromFile()
    pts1 = m1.points()

    m2 = modl.Model(model2)
    m2.loadFromFile()
    pts2 = m2.points()

    log.info("Comparison between two 3D model {0:s}: 1->{1:}pts 2->{2:}pts".format(label,len(pts1),len(pts2)))

    tr3D = rot3D

    res1 = numpy.Inf
    res2 = numpy.nan

    use_polynom = True
    use_polynom = use_polynom and len(pts1)>10 and len(pts2)>10

    use_res3d = True
    
    if len(pts1)>20 and len(pts2)>20:
        order = 3
    else:
        order = 2

    order = 3
#    order = 3
#    order = 4

    iter=0

    while (numpy.isnan(res2) or res2!=res1 or iter<miniter) and iter<maxiter: # Adjust the 3D transformation to accomodate the most points

        iter += 1

        index = comparePoints( pts1.copy(), pts2.copy(), tr3D )
        index1,index2 = index[:,0],index[:,1]

        nn = len(index1)

        if use_polynom and iter>=miniter:
            P = util.polynomial_regression(pts2[index2],pts1[index1],[order,order,order])
            tr3D = util.PolynomialTransformation3D(order,numpy.row_stack((P[0].coeff,P[1].coeff,P[2].coeff)))
        else: # Use affine transformation otherwise
            P = util.polynomial_regression(pts2[index2],pts1[index1],[1,1,1])
            tr3D = util.AffineTransformation3D(numpy.row_stack((P[0].coeff,P[1].coeff,P[2].coeff)))

        if use_res3d:
            res = numpy.sqrt(numpy.sum((pts1[index1]-tr3D.forward(pts2[index2]))**2,axis=1))
        else:
            res = numpy.sqrt(numpy.sum(((pts1[index1]-tr3D.forward(pts2[index2]))**2)[:,:2],axis=1))

        res1 = res2
        res2 = numpy.mean(res)

        res_ = numpy.where(res<dm)
        mean_ = numpy.mean(res_)
        std_ = numpy.std(res_)

        thrshold = min(dm,mean_+0.75*std_)

        index1 = index1[numpy.argwhere(res<thrshold)]
        index2 = index2[numpy.argwhere(res<thrshold)]

        log.info("Found {0:} of corresponding points. Keep {1:} (after {2:} iterations;threshold {3:.1f})".format(nn,index1.size,iter,thrshold))

#        # Second pass for calculating the transformation
#
#        res_ = res[numpy.argwhere(res<dm)]
#        index1_ = index1[numpy.argwhere(res<(numpy.mean(res_)+numpy.std(res_)))]
#        index2_ = index2[numpy.argwhere(res<(numpy.mean(res_)+numpy.std(res_)))]
#
#        pts1_ = numpy.squeeze(pts1[index1_])
#        pts2_ = numpy.squeeze(pts2[index2_])
#
#        print pts1_.shape
#
#        if use_polynom and iter>=miniter:
#            P = util.polynomial_regression(pts2_,pts1_,[order,order,order])
#            tr3D = util.PolynomialTransformation3D(order,numpy.row_stack((P[0].coeff,P[1].coeff,P[2].coeff)))
#        elif iter>=miniter: # Use affine transformation otherwise
#            P = util.polynomial_regression(pts2_,pts1_,[1,1,1])
#            tr3D = util.AffineTransformation3D(numpy.row_stack((P[0].coeff,P[1].coeff,P[2].coeff)))


    final_res = res[numpy.argwhere(res<dm)]
    log.info("Mean: {0:.1f}    Std:{1:.1f}".format(numpy.mean(final_res),numpy.std(final_res)))

#    for index,res_ in enumerate(res):
#        print "%i: %s" %(index,res_)
#
#    import sys
#    sys.exit(0)


    # Recalculate the relative affine transformation once the iterations are finished
    P = util.polynomial_regression(numpy.squeeze(pts2[index2]),numpy.squeeze(pts1[index1]),[1,1,1])
    tr3D = util.AffineTransformation3D(numpy.row_stack((P[0].coeff,P[1].coeff,P[2].coeff)))

    if work_directory!=None and label:

        if not os.path.lexists(work_directory): os.makedirs(work_directory)

        # Build the model

        model_path = os.path.join( work_directory, "{0:s}.modl".format(label) )
        m = modl.Model( model_path, template=m1 )

        o1 = m.addNewObject()
        c1 = o1.addNewContour()
        c1.points = pts1

        o2 = m.addNewObject()
        c2 = o2.addNewContour()
        c2.points = tr3D.forward(pts2)

        oc = m.addNewObject()   # Common points
        cc = oc.addNewContour()
        cc.points = numpy.squeeze(c1.points[index1] + c2.points[index2])/2.0

        m.save()

        if displayModel: os.system("imod -V {0:s}".format(output))

        # Build the residual image

        if use_polynom:
            res3d_png = os.path.join( work_directory, "res3d_{0:s}-order2.png".format(label) )
        else:
            res3d_png = os.path.join( work_directory, "res3d_{0:s}-order1.png".format(label) )

        pylab.figure()
        pylab.plot(res)
        pylab.xlabel(r"$\theta$")
        pylab.ylabel(r"Residual Distances between potential Markers")
        pylab.title(r"3D map %s" %label)
        pylab.savefig(res3d_png)

        if displayResidual: pylab.show()

    return index1, index2, tr3D


def compareMarkerModel( model1, slice1, model2, slice2, rotation, output="test.modl", doShow=True ):

    m1 = modl.Model(model1)
    m1.loadFromFile()

    markers1 = {}

    i1 = []
    pts1 = []

    for iobj1,obj1 in enumerate(m1.objects):
        pts1_ = numpy.row_stack([ct.points for ct in obj1.contours])
        index1 = numpy.where(pts1_[:,2]==slice1)
        pt1 = numpy.squeeze(pts1_[index1])
        log.info("obj: {0:}     {1:s}".format(iobj1,pt1))
        if len(pt1)!=0:
            i1.append(iobj1)
            markers1[iobj1] = pt1
            pts1.append(pt1)

    m2 = modl.Model(model2)
    m2.loadFromFile()

    markers2 = {}

    i2 = []
    pts2 = []


    for iobj2,obj2 in enumerate(m2.objects):
        pts2_ = numpy.row_stack([ct.points for ct in obj2.contours])
        index2 = numpy.where(pts2_[:,2]==slice1)
        pt2 = numpy.squeeze(pts2_[index2])
        log.info("obj: {0:}     {1:s}".format(iobj2,pt2))
        if len(pt2)!=0:
            i2.append(iobj2)
            markers2[iobj2] = pt2
            pts2.append(pt2)

#    pts1 = numpy.row_stack(markers1.values())
#    pts2 = numpy.row_stack(markers2.values())


    pts1 = numpy.row_stack(pts1)
    pts2 = numpy.row_stack(pts2)

    index = comparePoints( pts1, pts2, rotation )

    index1, index2 = index[:,0], index[:,1]

    for i1_,i2_ in zip(index1, index2):
        log.info("b:{0:}=c:{1:}".format(i1[i1_],i2[i2_]))

    pts1 = pts1[index1]
    pts2 = pts2[index2]

    # Save comparison into a model

    if output!=None:

        m = modl.Model( output, template=m1 )

        for pt1,pt2 in zip(pts1,pts2):
            o = m.addNewObject()
            c1 = o.addNewContour()
            c1.addPoint(*pt1)
            c2 = o.addNewContour()
            c2.addPoint(*pt2)

        m.save()

        m1_ = modl.Model(model1 + "_",template=m1)

        for index in index1:
            o = m1_.addObject(m1.objects[i1[index]])

        m1_.save()


        m2_ = modl.Model(model2 + "_",template=m2)

        for index in index2:
            o = m2_.addObject(m2.objects[i2[index]])

        m2_.save()


    if doShow and output!=None:

        os.system("imod -V %s" %output)

    return i1,i2


#    modfile1 = "/ncmirdata3/sphan/daniela/81780/txbr-align/bin4/grid2_sect1_G4K_dish7area4_TS2-V2-mSOG_exp7_7_11_1c.mod"
#    modfile2 = "/ncmirdata3/sphan/daniela/81780/txbr-align/bin4/grid2_sect1_G4K_dish7area4_TS2-V2-mSOG_exp7_7_11_1c.mod"
#
#    m1 = modl.Model(modfile1)
#    m1.loadFromFile()
#
#    XYZ1 = []
#
#    for index in index1:
#        obj1 = m1.objects[i1[index]]
#        XYZ1.append(obj1.points())
#
#    XYZ1 = numpy.row_stack(XYZ1)
#
#    m2 = modl.Model(modfile2)
#    m2.loadFromFile()
#
#    XYZ2 = []
#
#    for index in index2:
#        obj2 = m1.objects[i2[index]]
#        XYZ2.append(obj2.points())
#
#    XYZ2 = numpy.row_stack(XYZ2)
#
#    print XYZ1
#    print XYZ2




def comparePoints( pts1, pts2, transformation ):
    '''Compare two sets of (3d or 2d) points given a transformation.

    :param pts1: The first set of points
    :param pts2: The second set of points
    :param transformation: The transformation that roughly maps the two sets with
         *pts1=transformation(pts2)*
    '''

    in2D = False
    in3D = False

    if isinstance(transformation,util.PolynomialTransformation3D): in3D = True
    if isinstance(transformation,util.AffineTransformation2D): in2D = True

    if in2D:
        pts1_ = pts1[:,:2].copy()
        pts2_ = pts2[:,:2].copy()
        pts2_ = transformation.forward(pts2_)

    if in3D:
        pts1_ = pts1[:,:3].copy()
        pts2_ = pts2[:,:3].copy()
        pts2_ = transformation.forward(pts2_)

        pts1_ = pts1_[:,:2]
        pts2_ = pts2_[:,:2]

    # Try to match points between the two models

    tree1 = scipy.spatial.KDTree(pts1_)
    tree2 = scipy.spatial.KDTree(pts2_)
    
    d12,i12 = tree1.query(pts2_)
    d21,i21 = tree2.query(pts1_)

    check1 = (i12[i21]==range(len(pts1_)))
    check2 = (i21[i12]==range(len(pts2_)))

    index1 = numpy.where(check1==True)
    index2 = i21[index1]

    index1 = numpy.squeeze(index1)

    d = d21[index1]

    n = index1.size

    print "Found %i matching points!" %(n)
    print "Min/Max distance between two matching points: [%.1f,%.1f]" %(numpy.min(d),numpy.max(d))

    return numpy.column_stack( (index1, index2) )


def bckprjImageAtZ(  src, trf3D, Z=0.0 ):
    '''Backproject an image at Z=0 given the projection map is trf3D'''

    T = trf3D.T[:2] + trf3D.M[:2,2]*Z
    M = trf3D.M[:2,:2]

    warp2D = util.AffineTransformation2D( numpy.column_stack((T,M)) )
    warp2D = warp2D.inv()

    u = util.warp( src, warp2D )

    return numpy.where(numpy.isnan(u),0.0,u)


def scanProjectSeries( project, zoom=1.0, ncommon=None, dm=2.5, blur=False, reset=False  ):
    """This routine allows to find corresponding transformation between multiple series
    of the same area. It analyses the set of series within the project pair by pair,
    determine the corresponding 3D transformation between each of them. Prior to
    this routine, the series must have been aligned individually; the tracks, projection
    maps and 3D markers positions must be available as input files with extension '.indiv'
    (i.e. )

    :param project: The TxBR project made of multiple series.
    :param zoom: Zoom factor.
    :param ncommon: The procedure will only keep markers that have been detected in
        at lest ncommon series. If set to None, keep all the markers.
    :param reset: Re-evaluate the rotation angle between series if needed.

    """

    print "Scan the project %s" %(project)


    if len(project.series)<=1:
        print "The project does not contain multiple tilt series."
        return

  #  zoom=0.75

    restrict2common = ncommon!=None

    apply_gaussian_filter = blur
    gaussian_filter_ker_size = 2.0
    work_on_mid_tomogram = False

    counter = 0
    map = {}
    relativeTransforms =  __loadTransformations__(project)

    l0 = len(project.basename)

    # For each pair of series, try to determine the best affine transformation

    for index_ref,s_ref in enumerate(project.series):

        # HEEEERE
       # if project.extensions[index_ref]!='a': continue

        flip = False

        u_ref = s_ref.stack.getReferenceImage()
        if flip: u_ref = u_ref.T
        origin_ref = numpy.asarray(u_ref.shape)/2.0

        if work_on_mid_tomogram:
            trf3D_ref =  s_ref.projection.getAffineTransformation(s_ref.stack.getReferenceIndex())
            u_ref = bckprjImageAtZ(  u_ref, trf3D_ref )

        if apply_gaussian_filter: u_ref = scipy.ndimage.filters.gaussian_filter(u_ref,gaussian_filter_ker_size)

        mod_ref = os.path.join( project.align_directory, "%s.mod.indiv" %(s_ref.basename) )

        for index,s in enumerate(project.series): 

            if index<=index_ref: continue

            label = "%s-%s" %(s_ref.basename[l0:],s.basename[l0:]) # A label for the final plot images
            label_ = "%i-%i" %(index_ref,index)

            # Step (i): Compare the reference projection images to assert a relative 3D transformation
            # mapping the two series

            if reset or not label_ in relativeTransforms:

                u = s.stack.getReferenceImage()

                if work_on_mid_tomogram:
                    trf3D =  s.projection.getAffineTransformation(s.stack.getReferenceIndex())
                    u = bckprjImageAtZ(  u, trf3D )

                if apply_gaussian_filter: u = scipy.ndimage.filters.gaussian_filter(u,gaussian_filter_ker_size)

                rot = scanrot( u_ref, u, angles=(0,2*numpy.pi,150), zoom=zoom, label=label, work_directory=project.multser_dir, display_slices=False, display_peak=False )
   #             rot = scanrot( u_ref, u, angles=(0,2*numpy.pi,360), zoom=zoom, label=label, work_directory=project.multser_dir, display_slices=False, display_peak=False )

                rotinv = rot.inv()
                origin = rotinv.forward(origin_ref)

                angle = evaluateAngle( u_ref, origin_ref, u, origin )
                print "Matching angle from polar coordinates cross-correlation: %.1f" %(numpy.degrees(angle))

                t3D = numpy.append(rot.T,0.0)
                m3D = numpy.row_stack((rot.M,(0.0,0.0)))
                a = numpy.array(((0.0),(0.0),(1.0)))
                m3D = numpy.column_stack((m3D,a))

                rot3D = util.AffineTransformation3D(numpy.column_stack((t3D,m3D)))

            else:

                rot3D = relativeTransforms[label_]

#            print rot3D.T
#            print rot3D.M

            if flip:
                temp = rot3D.T[0]
                rot3D.T[0] = rot3D.T[1]
                rot3D.T[1] = temp

                rot3D.M[:2,:2] = rot3D.M[:2,:2].T

#            print rot3D.T
#            print rot3D.M

            
            # Step (ii): Find the corresponding markers between the series, and refine the relative
            # 3D transformation

            mod = os.path.join( project.align_directory, "%s.mod.indiv" %(s.basename) )

            index1, index2, tr3D = compare3DModel( mod_ref, mod, rot3D, dm=dm, label=label, work_directory=project.multser_dir, displayModel=False )

            relativeTransforms[label_] = tr3D
            
            for i1,i2 in zip(index1,index2):
                key1 = "%i-%i" %(index_ref,i1)
                key2 = "%i-%i" %(index,i2)
                if not key1 in map and not key2 in map:
                    map[key1] = counter
                    map[key2] = counter
                    counter += 1
                elif key1 in map and not key2 in map:
                    map[key2] = map[key1]
                elif not key1 in map and key2 in map:
                    map[key1] = map[key2]
                elif map[key1]!=map[key2]:
                    index1,index2=map[key1],map[key2]
                    for key,value in map.iteritems():
                        if value==index2: map[key] = index1

            try: # Save the relative transforms at each step
                __saveTransformations__( project, relativeTransforms )
            except:
                print "Unexpected error:", sys.exc_info()[0]
                pass

    # Reorder the constraints map

    markers = [[] for i in range(counter)]
    for key in map:
        markers[map[key]].append(key)

    pbs = []
    markers_ = []

    imarker = 0
    for index,m in enumerate(markers): # Merge the potential duplicate
        if len(m)==0: continue
        print "#%i ( obj %i - %i pts) %s" %(imarker,index,len(m),m)
        tilts = [item.split("-")[0] for item in m]
        counts = collections.Counter(tilts)
        if len(counts)!=len(tilts):
            pbs.append(m)
            continue
        imarker += 1
        markers_.append(m)

    markers = markers_

    print "There are %i common markers" %(imarker+1)

    if len(pbs)>0:
        print "Problem with %i clusters:" %(len(pbs))
        print pbs
    else:
        print "No apparent problem:"

    if restrict2common: # Restrict to markers present in every tilt series.
        markers = [ m for m in markers if len(m)>=ncommon]

    # Create the marker list and the 3D model

    markerfiles_in = []
    markerfiles_out = []
    models_3D_in = []
    XYZ = []

    for s in project.series:
        file_in = os.path.join( project.align_directory, "%s.mrk.indiv" %(s.basename) )
        m_in = modl.Model(file_in)
        m_in.loadFromFile()
        markerfiles_in.append(m_in)
        file_out = os.path.join( project.align_directory, "%s.mrk" %s.basename )
        m_out = modl.Model(file_out,template=m_in)
        markerfiles_out.append(m_out)
        mod3D_in = os.path.join( project.align_directory, "%s.mod.indiv" %(s.basename) )
        mod3D_in = modl.Model(mod3D_in)
        mod3D_in.loadFromFile()
        models_3D_in.append(mod3D_in)

    for index,marker in enumerate(markers):
        pts = []
        for iseries in range(len(project.series)):
            markerfiles_out[iseries].addNewObject()
        for token in marker:
            indexOfSeries,indexOfMarker = token.split("-")
            indexOfSeries = int(indexOfSeries)
            indexOfMarker = int(indexOfMarker)
            markerfiles_in[indexOfSeries]
            markerfiles_in[indexOfSeries].objects[indexOfMarker]
            markerfiles_out[indexOfSeries]
            markerfiles_out[indexOfSeries].objects[index]
            markerfiles_out[indexOfSeries].objects[index].contours = markerfiles_in[indexOfSeries].objects[indexOfMarker].contours
            pt3d = numpy.squeeze(models_3D_in[indexOfSeries].objects[indexOfMarker].points())
            if indexOfSeries!=0:
                tr3D = relativeTransforms["%i-%i" %(0,indexOfSeries)]
                pt3d = tr3D.forward(pt3d)
            pts.append(pt3d)
        XYZ.append(pts)

    output_model = os.path.join( project.multser_dir, "full.modl" )
    modl.saveAs3DModel( output_model, [XYZ] )

    if not restrict2common:
        for indexOfSeries,m in enumerate(markerfiles_in):
            n = 0
            for indexOfMarker,obj in enumerate(m.objects):
                key = "%i-%i" %(indexOfSeries,indexOfMarker)
                if key in map: continue # Not a singleton
                n += 1
                for iseries in range(len(project.series)):
                    markerfiles_out[iseries].addNewObject()
                markerfiles_out[indexOfSeries].objects[-1].contours = obj.contours
            print "%i aditional points for series %s" %(n,project.series[indexOfSeries].basename)

    for m in markerfiles_out:
        print "%s: %i points" %(m.filename,len(m.objects))
        m.save()

    # Correct the the projection maps so the different volumes match

    for iseries,s in enumerate(project.series):
        if iseries==0: continue
        tr3D = relativeTransforms["%i-%i" %(0,iseries)]
        tr3D = tr3D.inv()
        for index in range(s.numberOfExposures()):
            P = s.projection.getAffineTransformation(index)
            P = P.compose(tr3D)
            s.projection.x_coefficients[index,0] = P.T[0]
            s.projection.x_coefficients[index,1:4] = P.M[0,:]
            s.projection.y_coefficients[index,0] = P.T[1]
            s.projection.y_coefficients[index,1:4] = P.M[1,:]

#    print project.series[0].projection.x_coefficients
#    print project.series[0].projection.y_coefficients
#    print project.bin
#
#    import sys
#    sys.exit(0)

    project.saveAll()

    
if __name__ == '__main__':

    # To generate a seed file from one tilt to the other

    stfile1 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80182/processed_data/a/Ani12Gr3Syn2a.st"
    stfile2 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80182/processed_data/d/Ani12Gr3Syn2d.st"

    beadfile1 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80182/processed_data/a/Ani12Gr3Syn2a.seed"
    beadfile2 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80182/processed_data/d/Ani12Gr3Syn2d.seed"

    slice1 = 56
    slice2 = 18


    stfile1 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80092/processed_data/individuals/a/Ani12Gr2Syn1-CCDa.preali"
    stfile2 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80092/processed_data/individuals/c/Ani12Gr2Syn1-CCDc.preali"

    beadfile1 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80092/processed_data/individuals/a/Ani12Gr2Syn1-CCDa.seed"
    beadfile2 = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20162/microscopy_80092/processed_data/individuals/c/Ani12Gr2Syn1-CCDc.seed"

    slice1 = 56
    slice2 = 18


    # Rotate the second image and ross-correlate to find the matching angle

    f1 = mrc.MRCFile(stfile1)
    u1 = f1.getZSliceAt(slice1)

    f2 = mrc.MRCFile(stfile2)
    u2 = f2.getZSliceAt(slice2)

    # May be apply a gaussian filter
    
    u1 = scipy.ndimage.filters.gaussian_filter(u1,2)
    u2 = scipy.ndimage.filters.gaussian_filter(u2,2)

    # Find the best corresponding angle by cross-correlation

    n = float(500)
    zoom = n/float(max(f1.nx,f1.ny))

    rot = scanrot( u1, u2, angles=(0,2*3.14,300),zoom=zoom,display_slices=True )

    # Use polar coordinate images around a matchng point to get the matching angle

    t2d = rot.inv()

    origin1 = numpy.asarray(u1.shape)/2.0
    origin2 = t2d.forward(origin1)

    angle = evaluateAngle( u1, origin1, u2, origin2 )

    print "Matching angle from polar coordinates cross-correlation: %.1f" %(numpy.degrees(angle))

    rotateModel( beadfile1, t2d, output=beadfile2 )

 