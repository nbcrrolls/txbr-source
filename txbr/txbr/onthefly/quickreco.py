import sys
import os
import os.path
import shutil
import signal
import errno
import glob
import cPickle
import shutil
import multiprocessing
import numpy
import scipy.ndimage.interpolation
import scipy.stats
import mrc
import txbr.prefil
import txbr.utilities
import txbr.onthefly
import util
import modl
import phasecorr

global MESSAGE_FILE
global PID_FILE
global STATUS_FILE
global FILTER_DIR
global VOL_DIR
global CC_CROP_RATE
global Z_ASPECT_RATIO

MESSAGE_FILE = "message.txt"
PID_FILE = "pid.txt"
STATUS_FILE = "status.fly"

from txbr.txbrconf import DETECT_DIR
from txbr.txbrconf import ALIGN_DIR
from txbr.txbrconf import FILTER_DIR
from txbr.txbrconf import VOL_DIR
from txbr.txbrconf import WEB_DIR

from txbr import log

CC_CROP_RATE = 0.0
Z_ASPECT_RATIO = 0.1


def detect_task( stack, index, detect_dir, overwrite, nmin, nmax, size ):
    """
    The core task for detecting beads in images. Used with the multiprocessing module.
    """

    nx,ny = stack.getImageSize()

    detect_tif = os.path.join(detect_dir,"{0:s}.detect.{1:03}.tif".format(stack.basename,index))
    detect_mod = os.path.join(detect_dir,"{0:s}.bead.{1:03}.mod".format(stack.basename,index))
    detect_txt = os.path.join(detect_dir,"{0:s}.all.bead.{1:03}.txt".format(stack.basename,index))

    import scipy.ndimage.filters
    if not overwrite and os.path.exists(detect_tif) and os.path.exists(detect_mod):
        pts = modl.loadAllPoints(detect_mod)
    else:
        u = stack.getImageAt(index)

     #   u = scipy.ndimage.filters.gaussian_filter(u,1.5)

        pts, detector, patches = txbr.onthefly.findMarkersInImage( u, nmin=nmin, nmax=nmax, size=size, saveAllPointsInFile=detect_txt )
        pts = numpy.asarray(pts)
        Z = numpy.zeros((pts.shape[0],1))
        pts = numpy.append(pts,Z,axis=1)
        util.saveImage( detect_tif, detector, normalize=True )
        modl.saveAs3DModel( detect_mod, pts, nx, ny, 1, doShow=False )

    pts[:,2] = float(index)

    return pts


def detect_task_star( args ):
    """
    A wrapper to the *detect_task* function that takes a singe list argument.
    """

    return detect_task(*args)


def __bckprj_image_at__( src, trf3D, Z=0.0 ):
    '''
    Backproject an image at Z=0 given the projection map is trf3D.
    '''

    T = trf3D.T[:2] + trf3D.M[:2,2]*Z
    M = trf3D.M[:2,:2]

    warp2D = util.AffineTransformation2D( numpy.column_stack((T,M)) )
    warp2D = warp2D.inv()

    u = util.warp( src, warp2D )

    return numpy.where(numpy.isnan(u),0.0,u)


def __corr__( stack, angle, angle_ref, prj_id=None, correlate_before_stretch=True, gaussian_peak=True ):
    """
    Correlate two images from an image stack (between *angle* and *angle_ref*).

    :return: The affine transformation between the two corresponding images.
    """

    if angle==angle_ref:
        return util.Translation3D(0.0,0.0,0.0)

    if correlate_before_stretch:
        mat = stack.getImageAtAngle(angle)
        mat_ref = stack.getImageAtAngle(angle_ref)
    elif prj_id!=None:
        src = stack.getImageAtAngle(angle)
        src_ref = stack.getImageAtAngle(angle_ref)
        mat = __bckprj_image_at__( src, prj_id[index] )
        mat_ref = __bckprj_image_at__( src_ref, prj_id[index] )

    crop_rate = correlate_before_stretch and CC_CROP_RATE or max(0.1,CC_CROP_RATE)

    if gaussian_peak:
        tx,ty,widthx,widthy,theta,a,b =  phasecorr.correlate( mat_ref, mat, crop_rate=crop_rate )
        angleOfTrace.append([ widthy/widthx, theta ])
        log.info("Tilt sequence {0:3}/{1:-3}: t=({2:+7.3f},{3:+7.3f})  width=({4:+7.3f},{5:+7.3f}) ratio:{6:.3f}  angle={7:-3.3f}".format(index_ref,index,tx,ty,widthx,widthy,widthy/widthx,theta))
    else:
        tx,ty = util.crossCorrelate( mat_ref, mat, crop_rate=crop_rate)
        log.info("Tilt sequence {0:3}/{1:-3}: t=({2:+7.3f},{3:+7.3f})".format(index_ref,index,tx,ty))

    return util.Translation3D(tx,ty,0.0)


class QuickReconstruction():
    '''
    A class to perform very quick reconstruction on small volume, either once all
    the micrographs have been collected or on the fly.
    '''

    def __init__( self, stack, workDirectory=None, scope=None, cross_correlation_before_stretch=True,
                  START_FROM_SCRATCH = False ):

        self.stack = stack
        
        self.basename = stack.basename

        if workDirectory==None:
            self.workDirectory = self.stack.workDirectory
        else:
            self.workDirectory = workDirectory

        if scope==None:
            self.scope = self.stack.scope
        else:
            self.scope = scope

        self.showVolume = False
        self.makeAnimatedGif = True

        self.cross_correlation_before_stretch = cross_correlation_before_stretch
        self.START_FROM_SCRATCH = START_FROM_SCRATCH

        self.data_dir = os.path.join(self.stack.workDirectory,"data","bin{0:}".format(self.stack.bin))
        self.detect_dir = os.path.join( self.workDirectory, DETECT_DIR, "bin{0:}".format(self.stack.bin))
        self.align_dir = os.path.join( self.workDirectory, ALIGN_DIR, "bin{0:}".format(self.stack.bin))
        self.filt_dir = os.path.join( self.workDirectory, FILTER_DIR, "bin{0:}".format(self.stack.bin))
        self.vol_dir = os.path.join( self.workDirectory, VOL_DIR, "bin{0:}".format(self.stack.bin))
        #self.web_dir = os.path.join( self.workDirectory, WEB_DIR)
        self.web_dir = os.path.join( self.workDirectory, VOL_DIR)

        self.message_file = os.path.join( self.workDirectory, MESSAGE_FILE)
        self.pid_file = os.path.join( self.workDirectory, PID_FILE )
        self.status_file = os.path.join( self.workDirectory, STATUS_FILE)
        self.config_file = os.path.join( self.workDirectory, txbr.utilities.CONFIG_FILE)

        self.load()

        self.setRotationAxis( txbr.utilities.loadRotationAxis(self.scope) )


    def update_status(self):
        '''Update the status file'''

        f = open( self.message_file, 'wb' )
        f.write( self.status() )
        f.close();


    def status(self):
        '''Gives a status on the on the fly reconstruction'''

        try:
            try:
                f = open(self.pid_file, "rb")
                pid = int(f.read())
                f.close()
            except:
                raise OSError(errno.ESRCH, "No PID file")
            os.kill(pid, signal.SIG_DFL)
        except OSError, err:
            if err.errno == errno.ESRCH:
                return "Txbr on the fly is not running"

        return "Txbr on the fly is running on PID: {0:}".format(pid)
    

    def reset( self ):
        '''Kill the process and kill the files'''

        try:
            f = open(self.pid_file, "rb")
            os.kill(int(f.read()), signal.SIGKILL)
            f.close()
        except Exception, err:
            log.warning(err)

        self.clear()


    def clear(self):
        """
        Clear the file...
        """

        log.info("Clear...")

        files2remove = glob.glob(os.path.join(self.workDirectory,"*.fly"))
        files2remove.append(os.path.join(self.workDirectory,"result.html"))
        files2remove.append(os.path.join(self.workDirectory,"X-anim.gif"))
        files2remove.append(os.path.join(self.workDirectory,"Y-anim.gif"))
        files2remove.append(os.path.join(self.workDirectory,"Z-anim.gif"))

        for file in files2remove:
            if os.path.exists(file):
                os.remove(file)

        if os.path.exists(self.status_file):
            os.remove(self.status_file)

        if os.path.exists(self.message_file):
            os.remove(self.message_file)

        if os.path.exists(self.pid_file):
            os.remove(self.pid_file)

        if os.path.exists(self.filt_dir):
            shutil.rmtree(self.filt_dir)

        if self.vol_dir!="." and os.path.exists(self.vol_dir):
            shutil.rmtree(self.vol_dir)
        elif os.path.exists("{0:s}_vol.mrc".format(self.stack.basename)):
            os.remove("{0:s}_vol.mrc".format(self.stack.basename))

        self.load()


    def update( self, doBeadDetection=False, doReconstruction=True ):
        '''Update the reconstruction'''
        
        if doBeadDetection: self.updateBeadDetection()
        if doReconstruction: self.updateReconstruction()


    def updateBeadDetection( self, size=1000, factor=4 ):
        """
        Update the detection information
        """

        log.info("Update detection")

        try:
            (nx,ny) = self.stack.getImageSize()
            nx *= self.stack.bin
            ny *= self.stack.bin
            factor = max(1,int(round(max(nx,ny)/size)))
        except:
            pass

        files = self.stack.scanSourceDirectory( verbose=True )

        data_dir = os.path.join(self.stack.workDirectory,"data","bin{0:}".format(factor))
        detect_dir = os.path.join( self.workDirectory, DETECT_DIR, "bin{0:}".format(factor))

        if not os.path.lexists(data_dir):
            os.makedirs(data_dir)

        if not os.path.lexists(detect_dir):
            os.makedirs(detect_dir)

        for file in files:
            head,tail = os.path.split(file)
            basename,extension = os.path.splitext(tail)
            srcbin_tif = os.path.join(data_dir,"{0:s}.tif".format(basename))
            detect_tif = os.path.join(detect_dir,"{0:s}.detect.tif".format(basename))
            detect_mod = os.path.join(detect_dir,"{0:s}.bead.mod".format(basename))
            if os.path.lexists(srcbin_tif): continue
            open(srcbin_tif,'a').close()  # Create the new file
            child_pid = os.fork()
            if True and child_pid==0: # Child process
                log.info("Child Process: PID# {0:s}".format(os.getpid()))
                u = scipy.ndimage.interpolation.zoom(util.loadImage(file),1.0/factor)
                util.saveImage( srcbin_tif, u )
                points, detector, patches = txbr.onthefly.findMarkersInImage( u, nmax=100, size=11 )
                points = numpy.asarray(points)
                Z = numpy.zeros((points.shape[0],1))
                points = numpy.append(points,Z,axis=1)
                nx,ny = u.shape
                nz = 1
                util.saveImage( detect_tif, detector, normalize=True )
                modl.saveAs3DModel( detect_mod, points, nx, ny, nz, doShow=False )
                sys.exit(0) # Exit from the child process


    def finalizeBeadAlignement( self, size=1000, factor=4 ):

        try:
            (nx,ny) = self.stack.getImageSize()
            factor = int(round(max(nx,ny)/size))
        except:
            pass

        align_dir = self.align_dir + "_bin-{0:}".format(factor)
        pattern = os.path.join( align_dir, '*.bead.mod' )

        slices = glob.glob( pattern )
        slices.sort()

        points = []
        for slice in slices:
            pts = modl.loadAllPoints(slice)
            points.append(pts)

        points = []
        prjmap = []

    #    qa = txbr.onthefly.QuickAlign( basename, points, prjmap, xmin, xmax, ymin, ymax )


    def updateReconstruction(self):
        '''Update the reconstruction'''

        try:
            try:
                f = open(self.pid_file, "rb")
                pid = int(f.read())
                f.close()
            except:
                raise OSError(errno.ESRCH,"No PID file")
            os.kill(pid, 0) # Check is the process exists
        except OSError, err:
            if err.errno == errno.ESRCH:
                log.warning("Txbr on the fly is not running")
                f = open(self.pid_file, "wb")
                f.write('')
                child_pid = os.fork()
                if child_pid==0: # Child process
                    log.info("Child Process: PID# {0:s}".format(os.getpid()))
                    f.write("{0:}".format(os.getpid()))
                    f.close()
                    self.load()
                    while len(self.angles2process)!=0:
                        self.filter()
                        self.align()
                        self.bckprj()
                        self.save()
                        self.load()
                    os.remove(self.pid_file)
                    sys.exit(0) # Exit from the child process
            elif err.errno == errno.EPERM:
                log.warning("No permission to signal this process!")
            else:
                log.warning("Unknown error")
        else:
            log.warning("Txbr on the fly is already running")  # do nothing

        
    def load(self):
        '''
        Load the data
        '''

        self.stack.load()   # Scan for all images

        self.prjmap_id = {}
        self.prjmap = {}

        # Estimate the dimensions of the volume

        if not self.stack.isEmpty():

            self.nx, self.ny = self.stack.getImageSize()

            self.Zmax = Z_ASPECT_RATIO*max(self.nx,self.ny)
            self.Zmin = - self.Zmax
            self.nz = int(self.Zmax-self.Zmin)

        # What is to reconstruct?

        self.processed, self.processed_t = [],[]
        self.angles2process, self.index2process, self.tilts2process = [],[], []

        if ( not self.START_FROM_SCRATCH and os.path.exists(self.status_file) ):
            self.processed,self.processed_t = cPickle.load( open(self.status_file,'rb') )
        else:
            self.processed,self.processed_t = [],[]

        # Remove the data that is already processed

        self.angles2process = [ angle for angle in self.stack.angles if not angle in self.processed ]
        self.index2process = [ self.stack.indexOfAngle(angle) for angle in self.angles2process ]
        self.tilts2process = [ self.stack.slice4angle[angle] for angle in self.angles2process ]

        # Some display

        log.info(self.stack.angles)
        log.info("Angles to be processed: {0:s}".format(self.angles2process))
        log.info("Processed angles: {0:s}".format(self.processed))


    def getFrame( self ):
        '''
        Get the reconstruction frame...
        '''

        origin = [1.0,1.0,self.Zmin]
        end = [self.nx,self.ny,self.Zmax]

        return numpy.row_stack((origin,end))
    

    def saveAsTxBRProject( self ):
        '''
        Save this rough alignment in a txbr-configuration file.
        '''

        log.info("Save as a TxBR project")

        project = txbr.TxBRproject( self.stack.directory, self.stack.basename, self.stack.workDirectory, bin=self.stack.bin )

        project.reconstruction.setOrigin( 1.0, 1.0, int(self.Zmin) )
        project.reconstruction.setEnd( self.nx, self.ny, int(self.Zmax) )
        
        sx,sy = self.stack.getImageScale()
        project.reconstruction.scale.x = sx
        project.reconstruction.scale.y = sy
        project.reconstruction.scale.z = sx

        project.tiltAngles = self.stack.getTiltAngles()

        project.addSeries( extensionOfSeries='' )

        series = project.series[0]

        series.tiltAngles = self.stack.getTiltAngles()
        series.projection = txbr.ProjectionMap(1)

        series.sx = sx
        series.sy = sy

        ntilt = len(self.prjmap)

        series.prealignTranslations = numpy.zeros((ntilt,2))

        series.projection.x_coefficients = numpy.zeros((ntilt,4))
        series.projection.y_coefficients = numpy.zeros((ntilt,4))
        series.projection.scaling_coefficients = numpy.zeros((ntilt,4))

        for index,p in self.prjmap.iteritems():
            series.projection.x_coefficients[index,0] = p.T[0]
            series.projection.x_coefficients[index,1:] = p.M[0,:]
            series.projection.y_coefficients[index,0] = p.T[1]
            series.projection.y_coefficients[index,1:] = p.M[1,:]
            series.projection.scaling_coefficients[index,0] = 1

        series.rotAxis = self.rot_axis
        #project.reconstruction.rotAxis = self.rot_axis
        
        project.save()


    def setProjectionMap( self, prjmap ):
        '''
        Set the projection maps (in the case there are already know).
        Variable "maps" should be a dictionary containing the index of a view as 
        a key and the corresponding PolynomialTRansformation3D as a value
        '''

        self.prjmap = prjmap


    def setRotationAxis( self, rot_axis ):
        '''
        Set the rotation axis and the filter angle.
        '''

        self.rot_axis = rot_axis
        
        if self.rot_axis[1]==0:
            self.filt_angle = numpy.degrees(numpy.pi/2.0)
        else:
            self.filt_angle = numpy.degrees(numpy.arctan(-self.rot_axis[0]/self.rot_axis[1]))

        self.filt_angle = -self.filt_angle


    def align( self, gaussian_peak=True, refine_rot_axis=False, saveAsTxBRProject=False ):
        '''In this "align" routines, the images are aligned between themselves.
        If a projection correspondence is not available in the "prjmap" dictionary,
        then an image is cross-correlated with a reference frame. In this case, the
        images can be cross-correlated relatively to the first tilt or sequentially.
        '''
        
        if len(self.angles2process)==0: return

        cc_seq = True    # If True, images are cross-correlated sequentially
        from_scratch = len(self.processed)==0   # If True, no images have already been processed

        log.info("Alignment: sequential={0:}   from scratch={1:}".format(cc_seq,from_scratch))

        if not from_scratch:
            if cc_seq:
                index_ref = len(self.processed)-1
                angle_ref = self.processed[-1]
            else:
                index_ref = 0
                angle_ref = self.processed[0]
        else:
            if cc_seq:
                index_ref = 0
                angle_ref = self.angles2process[0]
            else:
                index_ref = self.stack.getReferenceIndex()
                angle_ref = self.stack.getReferenceAngle()

        # Define the ideal projection map

        self.prjmap_id = {}

        for index,angle in enumerate(self.angles2process):
            angle_in_rad = numpy.radians(angle)
            self.prjmap_id[index] = util.Rotation3D( (self.nx/2.0,self.ny/2.0,0.0), self.rot_axis, angle_in_rad )

        # Eventually backproject the reference image at Z=0

        if self.cross_correlation_before_stretch:
            matref = self.stack.getImageAt(index_ref)
        else:
            src = self.stack.getImageAt(index_ref)
            matref = __bckprj_image_at__( src, self.prjmap_id[index_ref] )

        # Calculate all the translations coefficients by cross-correlation

        prealignT = []
        angleOfTrace = []

        for index,angle in enumerate(self.angles2process):

            if angle==angle_ref:
                prealignT.append( util.Translation3D(0.0,0.0,0.0) )
                continue

            if self.cross_correlation_before_stretch:
                #mat = self.stack.getImageAt(index)
                mat = self.stack.getImageAtAngle(angle)
            else:
                #src = self.stack.getImageAt(index)
                src = self.stack.getImageAtAngle(angle)
                mat = __bckprj_image_at__( src, self.prjmap_id[index] )

            if self.cross_correlation_before_stretch:
                crop_rate = CC_CROP_RATE
            else:
                crop_rate = max(0.1,CC_CROP_RATE)

            #print "Crop rate: %f" %crop_rate

            if gaussian_peak:
                tx,ty,widthx,widthy,theta,a,b =  phasecorr.correlate( matref, mat, crop_rate=crop_rate )
                angleOfTrace.append([ widthy/widthx, theta ])
                log.info("Tilt sequence {0:3}/{1:-3}: t=({2:+7.3f},{3:+7.3f})  width=({4:+7.3f},{5:+7.3f}) ratio:{6:.3f}  angle={7:-3.3f}".format(index_ref,index,tx,ty,widthx,widthy,widthy/widthx,theta))
            else:
                tx,ty = util.crossCorrelate( matref, mat, crop_rate=crop_rate)
                log.info("Tilt sequence {0:3}/{1:-3}: t=({2:+7.3f},{3:+7.3f})".format(index_ref,index,tx,ty))

            prealignT.append( util.Translation3D(tx,ty,0.0) )

            if cc_seq:
                index_ref = index
                angle_ref = angle
                matref = mat

        t0x = 0.0
        t0y = 0.0
        
        if cc_seq:  # Transfer those parameter in the absolute frame
            if not from_scratch and len(self.processed_t)!=0:
                t0x = self.processed_t[-1][0]
                t0y = self.processed_t[-1][1]
        else:
            t0x = 0.0
            t0y = 0.0

        T0 = util.Translation3D(t0x,t0y,0.0)
        for index in range(len(prealignT)):
            prealignT[index] = prealignT[index].compose(T0)
            T0 = prealignT[index]

        for index in range(len(prealignT)):
            T = prealignT[index].inv()
            if self.cross_correlation_before_stretch:
                self.prjmap[index] = T.compose(self.prjmap_id[index])
            else:
                self.prjmap[index] = self.prjmap_id[index].compose(T)

        self.prealignT = prealignT

        # Evaluate the rotation axis from the traces

        angleOfTrace = numpy.asarray(angleOfTrace)

        if refine_rot_axis and angleOfTrace.size!=0:
            index = numpy.argsort(angleOfTrace[:,0],axis=0)
            ntrace = min(10,angleOfTrace.size)
            angleOfTrace = angleOfTrace[index]
            log.info("Angle of traces: {0:s}".format(str(angleOfTrace[-ntrace:,:])))
            #angleOfTrace_in_degrees = numpy.mean(angleOfTrace[-ntrace:,1])
            angleOfTrace_in_degrees = numpy.median(angleOfTrace[-ntrace:,1])
            angleOfTrace_in_radians = numpy.radians(angleOfTrace_in_degrees)
            axis = [ -numpy.sin(angleOfTrace_in_radians), numpy.cos(angleOfTrace_in_radians), 0.0 ]
            log.info("Calculated rotation axis: [{0:.3f},{1:.3f},{2:.3f}]; Angle: {3:.3f}".format(axis[0],axis[1],axis[2],angleOfTrace_in_degrees))
            self.setRotationAxis( axis )

        if saveAsTxBRProject:

            self.saveAsTxBRProject()
            
            txbr_file = os.path.join( self.stack.workDirectory, self.stack.basename + txbr.utilities.TXBR_EXTENSION )
            txbr_backup_file = os.path.join( self.stack.workDirectory, self.stack.basename + txbr.utilities.TXBR_EXTENSION + ".cross-corr.bin{0:}".format(self.stack.bin ))
            shutil.copy(txbr_file,txbr_backup_file)

        return self.prjmap


    def detect( self, size=6, nmin=0, nmax=10000, views=None, overwrite=False, doShow=False ):
        '''Detect the beads'''

        if not os.path.lexists(self.data_dir):
            os.makedirs(self.data_dir)

        if not os.path.lexists(self.detect_dir):
            os.makedirs(self.detect_dir)

        nx,ny = self.stack.getImageSize()

        if views==None:
            views = range(self.stack.getNumberOfImages())
            indexref = self.stack.getReferenceIndex()
        else:
            indexref = 0

        # First pass: Detect the markers

        pool = multiprocessing.Pool( processes=txbr.PROCESS_POOL_SIZE )

        tasks = [ [ self.stack, index, self.detect_dir, overwrite, nmin, nmax, size ] for index in views ]

        markers = pool.map( detect_task_star, tasks)

        pool.close()
        pool.terminate()

        log.info("indexref={0:}   npeak={1:}".format(views[indexref],len(markers[indexref])))

        # Second pass: homogeneize the number and global locations of the peaks
        # We are only interested in markers whick backproject in the volume

        P = util.AffineTransformation3D( numpy.array([[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0]]) )
        
        pref = P.compose(self.prjmap[indexref])

        BB = numpy.array([[0.0,0.0,0.0],[0.0,ny,0.0],[nx,0.0,0.0],[nx,ny,0.0]])

        for iview,m in enumerate(markers):  # Step 1. Determine the relevant number of peaks

            index = views[iview]

            p = P.compose(self.prjmap[index])
            p = P.compose(pref).compose(p.inv())

            pts = markers[iview]

            u = p.forward(pts)
            
            ii = numpy.where((u[:,0]>=0)*(u[:,0]<nx)*(u[:,1]>=0)*(u[:,1]<ny))
            ii = [range(len(u))]

            n = len(ii[0])

            if n==0:
                markers[iview] = numpy.zeros((0,3))
            else:
                markers[iview] = pts[ii[0]]
                markers[iview][:,2] = index

            if n!=0:
                min1 = numpy.min(u,axis=0)
                max1 = numpy.max(u,axis=0)
                u = p.forward(pts[ii[0]])
                min2 = numpy.min(u,axis=0)
                max2 = numpy.max(u,axis=0)
                log.debug("X: [{0:.1f},{1:.1f}] -> [{2:.1f},{3:.1f}] ".format(min1[0],max1[0],min2[0],max2[0]))
                log.debug("Y: [{0:.1f},{1:.1f}] -> [{2:.1f},{3:.1f}] ".format(min1[1],max1[1],min2[1],max2[1]))
                p =  P.compose(self.prjmap[index])
                p = P.compose(self.prjmap[index]).compose(p.inv())
                
            log.info("Number of interesting markers for tilt #{0:}: {1:}/{2:}".format(index,n,len(m)))

        peaks = numpy.array([ len(pts) for pts in markers ])
        npeak =  int(scipy.stats.cmedian(scipy.stats.trimboth(peaks,0.15)))
        npeak =  int(numpy.max(scipy.stats.trimboth(peaks,0.15)))
        log.info("npeak={0:} [Median value of the markers versus tilt]".format(npeak))

        for iview,m in enumerate(markers):  # Step 2. Select the appropriate number of peaks for each tilt

            index = views[iview]

            p =  P.compose(self.prjmap[index])
            p = pref.compose(p.inv())

            pts = markers[iview]

            u = p.forward(pts)
            ii = numpy.where((u[:,0]>=0)*(u[:,0]<nx)*(u[:,1]>=0)*(u[:,1]<ny))
            n = len(ii[0])

            if n<npeak:
                detect_txt = os.path.join(self.detect_dir,"{0:s}.all.bead.{1:03}.txt".format(self.stack.basename,index))
                pts = numpy.loadtxt(detect_txt)
                pts = numpy.asarray(pts)
                Z = numpy.ones((pts.shape[0],1))*float(index)
                pts = numpy.append(pts,Z,axis=1)
                u = p.forward(pts)
                ii = numpy.where((u[:,0]>=0)*(u[:,0]<nx)*(u[:,1]>=0)*(u[:,1]<ny))
                ii = [range(len(u))]
                n = len(ii[0])
                log.debug("Load more points [from pool of {0:}] for tilt {1:} -> {2:} in the right domain".format(len(pts),index,n))

            if n!=0:
                markers[iview] = pts[ii[0][:npeak]]
                markers[iview][:,2] = index
            else:
                markers[iview] = numpy.zeros((0,3))

            log.info("Number of markers for tilt #{0:}: {1:}/{2:}".format(index,n,len(markers[iview])))

        final_detect_mod = os.path.join(self.detect_dir,"{0:s}.bead.final.mod".format(self.stack.basename))
        modl.saveAs3DModel( final_detect_mod , markers, nx, ny, len(markers), doShow=doShow, mrcFileName=self.stack.mrc_file.filename )

        return markers
    

    def createStack( self, type="preali", Z=0, showVol=True ):
        '''Create a cross-correlated stack. There are two possibilities depending on
        the value of the variable *type*. If *type* is specified as "preali", then the
        correlated stack corresponding to an aligned tilt series. Otherwise, it corresponds
        to a stack where the slices would identify to the zero-tilt middle slice in the
        original volume (with no sharpening filter)'''

        if not os.path.lexists(self.align_dir):
            os.makedirs(self.align_dir)

        filename = os.path.join(self.align_dir,"{0:s}.preali".format(self.stack.basename))

        nx,ny = self.stack.getImageSize()
        nz = self.stack.getNumberOfImages()

        if len(self.prjmap)!=nz or len(self.prjmap_id)!=nz:
            self.align()

        log.info("Create the stack {0:} ({1:},{2:},{3:})".format(filename,nx,ny,nz))

        index_ref = self.stack.getReferenceIndex()

        if type=="preali":
            tf_ref = self.prjmap[index_ref].compose(self.prjmap_id[index_ref].inv())
        else:
            tf_ref = self.prjmap[index_ref]

        fout = mrc.MRCFile(filename)
        fout.setHeader(nx,ny,nz)

        for index in range(nz):

            src = self.stack.getImageAt(index)
            
            if type=="preali":
                tf = self.prjmap[index].compose(self.prjmap_id[index].inv())
                tf = tf.compose(tf_ref.inv())
            else:
                tf = self.prjmap[index]
                tf = tf.compose(tf_ref.inv())

            dest = __bckprj_image_at__( src, tf, Z=Z )

            fout.setZSliceAt(index,dest)

        fout.updateHeader()
        
        if showVol: os.system( "imod {0:s}".format(filename) )


    def filter(self):
        '''Filter the data'''

        log.info("Filter with angle: {0:.1f}".format(self.filt_angle))
        
        self.stack.filter( self.filt_angle, index2filter=self.index2process )


    def createVolume( self, output=None, from_scratch=True ):
        """Create the volume"""

        if not os.path.lexists(self.vol_dir):
            os.makedirs(self.vol_dir)

        log.info("Initialize MRC volume of size ({0:},{1:},{2:})".format(self.nx, self.ny, self.nz))

        if output==None:
            file = os.path.join( self.vol_dir, '{0:s}_vol.mrc'.format(self.stack.basename) )
        else:
            file = output

        mrcVolume = mrc.MRCFile(file)

        if from_scratch:
            mrcVolume.setHeader( self.nx, self.ny, self.nz)
            zeros = numpy.zeros((self.nx,self.ny))
            for iz in range(self.nz):
                mrcVolume.setZSliceAt(iz,zeros)

        return mrcVolume


    def bckprj( self, showVol=False ):
        '''Back-project the slices into the final volume'''

        if len(self.angles2process)==0: return

        mrcVolume = self.createVolume( from_scratch = len(self.processed)==0 )

        volumes = []

        for iz in range(self.nz):
            slice = mrcVolume.getZSliceAt(iz)
            volumes.append(slice)

        src_file = "{0:s}.SL".format(self.stack.basename)
        fsrc = mrc.MRCFile(src_file)

        for index,angle in enumerate(self.angles2process):

            sys.stdout.write( "Angle: {0:}: Z-> ".format(angle) )

        #    u = fsrc.getZSliceAt(index)
         #   u = self.stack.getFilteredImageAt(index)
            u = self.stack.getFilteredImageAtAngle(angle)

            for iz in range( self.nz ):

                Z = self.Zmin + iz

                mat = __bckprj_image_at__( u, self.prjmap[index], Z=Z )

                sys.stdout.write( '{0:}..'.format(Z) )
                sys.stdout.flush()

                volumes[iz] = mat + volumes[iz]

            sys.stdout.write("\n")

        for iz in range( self.nz ):
            u = numpy.asarray(volumes[iz])
            mrcVolume.setZSliceAt( iz, u )

        mrcVolume.updateHeader()

        if showVol: os.system("imod {0:s}_vol.mrc".format(self.stack.basename))


    def save(self):
        '''Save the data configuration'''

        if len(self.angles2process)==0: return

        # All the angles have been processed
        for angle,t in zip(self.angles2process,self.prealignT):
            if not angle in self.processed:
                self.processed.append(angle)
                self.processed_t.append(t.T[:2].tolist())

        f = open(self.config_file,'wb')
        f.write("{0:s}: [{1:.6f},{2:.6f},{3:.6f}]".format("Rotation Axis",self.rot_axis[0],self.rot_axis[1],self.rot_axis[2]))

        cPickle.dump((self.processed,self.processed_t),open(self.status_file,'wb'))

        mrcVolume = os.path.join( self.vol_dir, '{0:s}_vol.mrc'.format(self.stack.basename) )

        if os.path.lexists(mrcVolume) and self.showVolume:
            os.system( 'imod {0:s}'.format(mrcVolume) )

        if os.path.lexists(mrcVolume) and  self.makeAnimatedGif:
            if not os.path.lexists(self.web_dir):
                os.makedirs(self.web_dir)
            txbr.utilities.make_animated_gif( mrcVolume, self.web_dir )


if __name__ == '__main__':


    def onTheFlyAlignment( directory= "/home/sphan/data/MMV-7", basename = "MMV-7", cross_correlation_before_stretch=False, scope="j3200.de12" ):

        filename = "%s.st" %basename
        if not os.path.exists(filename):
            log.warning("File %s does not exist!" %filename)
            return

        stack = txbr.utilities.MRCStack("%s.st" %basename)
        reco = QuickReconstruction( stack, scope=scope, cross_correlation_before_stretch=cross_correlation_before_stretch )
        reco.align()
        reco.createStack( type="preali", showVol=True )
        reco.saveAsTxBRProject()

    def onTheFlyReconstruction( directory= "/home/sphan/data/MMV-7", basename = "MMV-7", cross_correlation_before_stretch=False, scope="j3200.de12" ):

        stack = txbr.utilities.MRCStack("{0:s}.st".format(basename))
        reco = QuickReconstruction( stack, scope=scope, cross_correlation_before_stretch=cross_correlation_before_stretch )
        reco.align()
        reco.filter()
        reco.bckprj( showVol=True)

    def directReconstruction( directory= "/home/sphan/data/MMV-7", basename = "MMV-7", scope="j3200.de12" ):

        import txbr

        project = txbr.TxBRproject( directory, basename )
        project.load()

        projection = project.series[0].projection

        prjmap = {}

        for index in range(projection.numberOfExposures()):
            prjmap[index] = projection.getAffineTransformation(index)


        stack = txbr.utilities.MRCStack("{0:s}.st".format(basename))
        reco = QuickReconstruction( stack, scope=scope )
        reco.setRotationAxis( project.series[0].rotAxis )
        reco.setProjectionMap( prjmap )
        reco.filter()
        reco.bckprj( showVol=True)

    # Different options

    options = { 1:onTheFlyReconstruction, 2:directReconstruction, 3:onTheFlyAlignment }

#    options[1]()
#    options[1]( basename = "FHV102709-29" )
#    options[2]( basename = "fhv6a" )
#    options[1]( basename = "MMV-9" )
#    options[1]( basename = "eel-syn_zero-loss_8k" )
#    options[1]( basename = "MMV-9" )
#    options[1]( directory="/ncmirdata5/sphan/temp/58930/bin8",basename = "FHV-3" )
#    options[1]( directory="/ncmirdata3/mterada/for_sebs",basename = "2Dbin8" )
#    options[1](directory="/ncmirdata5/sphan/temp/60111/bin8", basename="FHV-6")
#    options[3](directory="/ncmirdata3/sphan/otf", basename="test")
#    options[3](directory="/ncmirdata3/sphan/otf/14otf", basename="test")
#    options[1](directory="/ncmirdata3/sphan/otf/14otftiezt", basename="tietz", cross_correlation_before_stretch=False)
#    options[1](directory="/ncmirdata3/sphan/otf/15de12", basename="FHV-16_de12", scope="fei_titan.de12")
#    options[3](directory="/ncmirdata3/sphan/otf/15gatan", basename="FHV-16_gatan", scope="fei_titan")
#    options[3](directory="/ncmirdata3/sphan/FHV-18/bin4", basename="FHV-18a", scope="fei_titan.de12")
#    options[2](directory="/ncmirdata3/sphan/otf/mon", basename="SCN_night_OTO_cellchain001a", scope="j4000_2")
#    options[3](directory="/ncmirdata3/sphan/otf/mon2", basename="SCN_night_OTO_cellchain001a", scope="j4000_2")
#    options[3](directory="/ncmirdata3/sphan/otf/25", basename="FHV-25-2a", scope="fei_titan")


        
    def test( directory="/ncmirdata3/sphan/otf/ju", basename = "Ju_cont2.mrc", extension=txbr.utilities.TIFF, cross_correlation_before_stretch=False, scope="fei.titan" ):

        stack = txbr.utilities.ImageStack( basename, directory, extension=extension, scope="fei_titan")
#        stack.bin(2)

        reco = QuickReconstruction( stack, scope=scope, cross_correlation_before_stretch=cross_correlation_before_stretch,START_FROM_SCRATCH = True )
        reco.align()
      #  reco.createStack( type="preali", showVol=True )
        reco.filter()
        reco.bckprj( showVol=True)
       # reco.update()

    test()
 #   test(directory="/ncmirdata3/sphan/otf/ju2", basename = "CCDBid_74225_Ju_cont2.mrc", extension=txbr.utilities.JPG)
 