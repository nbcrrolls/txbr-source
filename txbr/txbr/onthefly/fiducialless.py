import sys
import os
import os.path
import signal
import errno
import glob
import re
import cPickle
import shutil
import numpy
import cv
import mrc
import txbr.prefil
import txbr.utilities

# 4000#1: Rotation Axis = [ -0.993129, 0.117027, 0.000000 ]

global SRC_DIR,FILT_DIR,VOL_DIR
global CONFIG_FILE,STATUS_FILE,PID_FILE
global E_X,E_Y,E_Z
global CC_CROP_RATE,Z_ASPECT_RATIO
global TIFF,JPG

SRC_DIR = "src_jpegs"
FILT_DIR = "filt_dir"
VOL_DIR = "vol_dir"
CONFIG_FILE= "config.fly"
MESSAGE_FILE = "message.txt"
PID_FILE = "pid.txt"
STATUS_FILE = "status.fly"

CC_CROP_RATE = 0.15
CC_CROP_RATE = 0.
Z_ASPECT_RATIO = 0.08

TIF = "tif"
JPG = "jpg"

E_X = ( 1.0, 0.0, 0.0 )
E_Y = ( 0.0, 1.0, 0.0 )
E_Z = ( 0.0, 0.0, 1.0 )

rotation_axis = "Rotation Axis"

def show( image ):
    '''Show an image with the openCV package'''
    
    cv.NamedWindow("debug")
    cv.ShowImage("debug", image)
    cv.WaitKey()
    

def img2mrc( file_in, file_out ):
    
    fout = mrc.MRCFile(file_out)
    
    src = cv.LoadImageM( file_in )
    ny,nx = cv.GetDims( src )
    dest = numpy.asarray( src, dtype='float' )
    dest =  dest[::-1,:,0].T
    
    fout.setHeader( nx, ny, 1 )
    fout.setZSliceAt( 0, dest )
    fout.updateHeader()
    
    return fout,(nx,ny)
    

def reset( basename, directory=SRC_DIR, extension=JPG, start_angle=-60.0, inc_angle=2.0 ):
    '''Reset the names of the tilt series'''

    print "Reset names for [basename] %s, [directory] %s, [extension] %s " %( basename,directory,extension )

    tilts = glob.glob( '%s/%s*.%s' %( directory, basename, extension ) )
    tilts.sort()

    print tilts
    
    for tilt in tilts:
        os.rename( tilt, '%s/%s_%.1f.%s' %( directory, basename, start_angle, extension ) )
        start_angle += inc_angle
    

class ReconstructionOnTheFly:
    '''Do a quick reconstruction on the fly...'''
    
    START_FROM_SCRATCH = False
    
    def __init__( self, directory, workDirectory, basename, extension=JPG, scope=None ):

        print "Reconstruction on the Fly"
        print "directory: %s" %directory
        print "workDirectory: %s" %workDirectory
        print "basename: %s" %basename
        print "scope: %s" %scope
        
        self.directory = directory
        self.workDirectory = workDirectory
        self.basename = basename
        self.extension = extension
        
        self.showVolume = False
        self.makeAnimatedGif = True

        if scope==None:
            config_file = CONFIG_FILE
        else:
            config_file = "config.%s.fly" %(scope.lower())

        data_dir =  os.path.join(os.path.dirname(__file__),'..','..','..','data')
        org_config_file = os.path.join(data_dir,'config',config_file)
        
        self.config_file = os.path.join( self.workDirectory, config_file)
        self.message_file = os.path.join( self.workDirectory, MESSAGE_FILE)
        self.pid_file = os.path.join( self.workDirectory, PID_FILE )
        self.status_file = os.path.join( self.workDirectory, STATUS_FILE)

        print "Copy %s to %s" % (org_config_file, self.config_file)
        shutil.copy(org_config_file, self.config_file)
        
        self.src_dir = os.path.join( self.directory, SRC_DIR)
        self.filt_dir = os.path.join( self.workDirectory, FILT_DIR )
        self.vol_dir = os.path.join( self.workDirectory, VOL_DIR )

        if not os.path.exists(self.workDirectory):
            os.mkdir(self.workDirectory)
        
        self.rot_axis = numpy.asarray(E_X)
        if self.rot_axis[1]==0:
            self.filt_angle = numpy.degrees(numpy.pi/2.0)
        else:
            self.filt_angle = numpy.degrees(numpy.arctan(-self.rot_axis[0]/self.rot_axis[1]))
        
        self.nx = 0
        self.ny = 0
        
        self.Zmax = Z_ASPECT_RATIO*numpy.max((self.nx,self.ny))
        self.Zmin = -self.Zmax
        self.nz = int(self.Zmax-self.Zmin)
        
        self.tilts2process = []


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
            
        return "Txbr on the fly is running on PID: %i" %pid
        

    def reset( self ):
        '''Kill the process and kill the files'''

        try:
            f = open(self.pid_file, "rb")
            os.kill(int(f.read()), signal.SIGKILL)
            f.close()
        except Exception, err:
            print err

        self.clear()


    def clear(self):

        files2remove = glob.glob(os.path.join(self.workDirectory,"*.fly"))

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
            
        if os.path.exists(self.vol_dir):
            shutil.rmtree(self.vol_dir)
        
        self.load()
        
        
    def update(self):
        '''Update the reconstruction'''

        try:
            try:
                f = open(self.pid_file, "rb")
                pid = int(f.read())
                f.close()
            except:
                raise OSError(errno.ESRCH,"No PID file")
            #print pid
            os.kill(pid, 0)
        except OSError, err:
            #print err
            if err.errno == errno.ESRCH:
                print "Txbr on the fly is not running"
                child_pid = os.fork()
                if child_pid == 0: # Child process
                    print "Child Process: PID# %s" % os.getpid()
                    f = open(self.pid_file, "wb")
                    f.write("%i" %os.getpid())
                    f.close()
                    self.load()
                    self.filter()
                    self.align()
                    self.bckprj()
                    self.save()
                    os.remove(self.pid_file)
            elif err.errno == errno.EPERM:
                print "No permission to signal this process!"
            else:
                print "Unknown error"
        else:
            print "Txbr on the fly is already running"  # do nothing

        
    def load(self):
        '''Load the data'''
        
        # The names of the file to process is a concatenation of a basename and the tilt angle
        
        pattern = os.path.join( self.src_dir, '%s_*.%s' %(self.basename,self.extension))
        print "Search for files like %s" %pattern
        self.tilts2process = glob.glob( pattern )
        self.unknownTilts = []
        
        # Extract the tilt angles
        
        self.angles2process = []
        
        for tilt in self.tilts2process:
            angle_match = re.match('\S*_([+-.\d]*).%s' %self.extension, tilt)
            if angle_match:
                angle = angle_match.group(1)
                self.angles2process.append( float(angle) )
            else:
                self.unknownTilts.append(tilt)
                
        self.angles2process.sort()
        
        self.tilts2process = [ tilt for tilt in self.tilts2process if not tilt in self.unknownTilts ]
        
        # Load some configuration information related to the used microscope, such as the
        # rotation axis. It will be used when filtering is applied.
        
        try:

            f = open(self.config_file, 'r')
            for line in f:
                rotAxis_match = re.match('%s:\s*\\[(\S*)\\]' %rotation_axis, line)
                if rotAxis_match:
                    rotAxis = rotAxis_match.group(1).split(',')
                    self.rot_axis = [float(tk) for tk in rotAxis]
            f.close()
            
        except IOError:
            
            print "Error loading the Configuration file!"
            
        if self.rot_axis[1]==0:
            self.filt_angle = numpy.degrees(numpy.pi/2.0)
        else:
            self.filt_angle = numpy.degrees(numpy.arctan(-self.rot_axis[0]/self.rot_axis[1]))
                    
        # Load some status information on the on the fly processing (How much of
        # the reconstruction is done...)
        
        if ( not self.START_FROM_SCRATCH and os.path.exists(self.status_file) ):
            self.processed,self.processed_t = cPickle.load( open(self.status_file,'rb') )
        else:
            self.processed,self.processed_t = [],[]
            
        # Remove the already processed data
        
        self.angles2process = [ angle for angle in self.angles2process if not angle in self.processed ]
        self.tilts2process = [ os.path.join(self.src_dir,'%s_%.1f.%s' %(self.basename,angle,self.extension)) for angle in self.angles2process ]
            
        # Some display
            
        #print "Files to be processed: %s" %(self.tilts2process)
        print "Unknown tilts: %s" %(self.unknownTilts)
        print "Angles to be processed: %s" %(self.angles2process)
        print "Processed angles: %s" %(self.processed)
        print "Rotation axis: %s" %(self.rot_axis)
        print "Filtering angle: %f" %(self.filt_angle)
        
        
    def filter(self):
        '''Filter the data'''
        
        if len(self.angles2process)==0: return
        
        filter_magnification = 1
        
        if not os.path.exists(self.filt_dir):
            os.mkdir(self.filt_dir)
            
        for file_in,angle in zip( self.tilts2process,self.angles2process ):
            file_mrc = "%s/%s_%.1f.mrc" %( self.filt_dir,self.basename,angle )
            file_out = "%s/%s_%.1f_f.mrc" %( self.filt_dir,self.basename,angle )
            print "Create file mrc %s..." %file_mrc
            f,(self.nx,self.ny) = img2mrc( file_in, file_mrc )  # First transfer the jpeg file into a mrc file
            self.Zmax = Z_ASPECT_RATIO*numpy.max((self.nx,self.ny))
            self.Zmin = -self.Zmax
            self.nz = int(self.Zmax-self.Zmin)
            txbr.prefil.rot_filter( file_mrc, file_out, filter_magnification, self.filt_angle)
        
        
    def align(self):
        '''In this function, images are cross-correlated.
        Images can be cross-correlated relatively to the first tilt or sequentially.
        ''' 
        
        if len(self.angles2process)==0: return
        
        cc_seq = True    # If True cross-correlate sequentially
                
        from_scratch = len(self.processed)==0
        
        if not from_scratch:
            if cc_seq:
                angle_ref = self.processed[-1]
            else:
                angle_ref = self.processed[0]
        else:
            angle_ref = self.angles2process[0]
        
        image_ref = os.path.join(self.src_dir,'%s_%.1f.%s' %( self.basename, angle_ref, self.extension ))
        image_ref =  cv.LoadImageM(image_ref)    
        matref = self.stretch_image( image_ref, numpy.radians(angle_ref) )
        
        self.nx = matref.cols
        self.ny = matref.rows
        
        self.Zmax = Z_ASPECT_RATIO*numpy.max((self.nx,self.ny))
        self.Zmin = -self.Zmax
        self.nz = int(self.Zmax-self.Zmin)
        
        # Calculate all the translations coefficients by cross-correlation
        
        t = numpy.zeros((len(self.angles2process),2))
            
        for index,angle in enumerate(self.angles2process):
            
            if from_scratch and angle==angle_ref: continue
        
            angle_in_rad = numpy.radians(angle)
            
            image = os.path.join(self.src_dir,'%s_%.1f.%s' %( self.basename, angle, self.extension ))
            mat = self.stretch_image( cv.LoadImageM( image ), angle_in_rad )
            
            t_ = self.cross_correlate( matref, mat )
            
            t[index,0] = t_[0]
            t[index,1] = t_[1]
            
            print "tile=%.1f (tx,ty)=(%.2f,%.2f)" %(angle,t[index,0],t[index,1])
            
            if cc_seq:
                angle_ref = angle
                matref = mat
          
        if cc_seq:  # Transfer thos parameter in the absolute frame
            t = numpy.cumsum(t,axis=0)
            if not from_scratch and len(self.processed_t)!=0:
                t[:,0] += self.processed_t[-1][0]
                t[:,1] += self.processed_t[-1][1]
              
        self.t = t
        
            
    def stretch_image( self, src, angle, Z=0.0, extra_translation=(0.0,0.0) ):
        '''Helper function to stretch an image (input data: MAT array) that was taken at a given tilt "angle" 
        for it to correspond to the slice at Z=0 in the final volume''' 
    
        axis = self.rot_axis
    
        nx,ny = cv.GetDims(src)
        
        dest = cv.CreateMat(nx,ny,cv.GetElemType(src))
        
        warp = cv.CreateMat(2,3,cv.CV_32FC1)
        
        Origin = cv.CreateMat(3,1,cv.CV_32FC1)
        u = cv.CreateMat(3,1,cv.CV_32FC1)
        
        Origin[0,0], Origin[1,0], Origin[2,0] = nx/2.0, ny/2.0, 0.0
        u[0,0], u[1,0], u[2,0] = axis[0], axis[1], axis[2]
        
        T = cv.CreateMat(2,1,cv.CV_32FC1)
        R = cv.CreateMat(2,2,cv.CV_32FC1)
    
        cos_phi = numpy.cos(angle)
        sin_phi = numpy.sin(angle)
    
        R[0,0] = u[0,0]*u[0,0]*(1-cos_phi) + cos_phi
        R[0,1] = u[0,0]*u[1,0]*(1-cos_phi) - u[2,0]*sin_phi
        
        R[1,0] = u[1,0]*u[0,0]*(1-cos_phi) + u[2,0]*sin_phi
        R[1,1] = u[1,0]*u[1,0]*(1-cos_phi) + cos_phi
        
        T[0,0] = Origin[0,0] - ( R[0,0]*Origin[0,0] + R[0,1]*Origin[1,0] )
        T[1,0] = Origin[1,0] - ( R[1,0]*Origin[0,0] + R[1,1]*Origin[1,0] )
                    
        T[0,0] -= ( u[0,0]*u[2,0]*(1-cos_phi) + u[1,0]*sin_phi )*Z
        T[1,0] -= ( u[1,0]*u[2,0]*(1-cos_phi) - u[0,0]*sin_phi )*Z
        
        # Take the inverse transformation
        
        Tinv = cv.CreateMat(2,1,cv.CV_32FC1)
        Rinv = cv.CreateMat(2,2,cv.CV_32FC1)
        
        cv.Invert( R, Rinv )
        Tinv[0,0] = - ( Rinv[0,0]*T[0,0] + Rinv[0,1]*T[1,0] )
        Tinv[1,0] = - ( Rinv[1,0]*T[0,0] + Rinv[1,1]*T[1,0] )
        
        # Apply the warp onto the image
        
        warp[0,0] = Rinv[0,0]
        warp[0,1] = -Rinv[0,1]
        warp[1,0] = -Rinv[1,0]
        warp[1,1] = Rinv[1,1]
        
        warp[0,2] = Tinv[0,0] + Rinv[0,1]*ny
        warp[1,2] = ny - Tinv[1,0] - Rinv[1,1]*ny
        
        warp[0,2] += extra_translation[0]
        warp[1,2] += extra_translation[1]
        
        cv.WarpAffine(src, dest, warp)
        
        return dest
    
    
    def cross_correlate( self, mat1, mat2 ):
        '''Cross-correlate two openCV MAT objects'''
        
        row_margin = int(CC_CROP_RATE*mat2.rows)
        column_margin = int(CC_CROP_RATE*mat2.cols)
    
       # mat2_ = cv.GetSubRect(mat2,(row_margin, column_margin, mat2.rows-2*row_margin, mat2.cols-2*column_margin))
        mat2_ = cv.GetSubRect(mat2,( column_margin, row_margin, mat2.cols-2*column_margin, mat2.rows-2*row_margin))
        template = cv.CreateMat(2*row_margin+1, 2*row_margin+1, cv.CV_32FC1)
        #template = cv.CreateMat(2*row_margin+1, 2*row_margin+1, cv.CV_64FC1)
        
        cv.MatchTemplate( mat1, mat2_, template, cv.CV_TM_CCORR_NORMED )
        
        ( minVal, maxVal, minLoc, maxLoc ) = cv.MinMaxLoc(template)
        t = numpy.array(maxLoc) - numpy.array((row_margin, column_margin))
        
        return t
    
    
    def createVolume( self, from_scratch=True ):
        
        if not os.path.lexists(self.vol_dir):
            os.mkdir(self.vol_dir)
            
        print "Initialize MRC volume of size (%i,%i,%i)" %(self.nx, self.ny, self.nz)
        
        file = os.path.join( self.vol_dir, '%s_vol.mrc' %( self.basename ) )
        
        mrcVolume = mrc.MRCFile(file)
        
        if from_scratch:
            mrcVolume.setHeader( self.nx, self.ny, self.nz)
            zeros = numpy.zeros((self.nx,self.ny))
            for iz in range(self.nz):
                mrcVolume.setZSliceAt(iz,zeros)
            
        return mrcVolume
    
    
    def bckprj(self):
        '''Back-project the slices into the final volume'''
        
        if len(self.angles2process)==0: return
        
        mrcVolume = self.createVolume( from_scratch = len(self.processed)==0 )
        
        # Process the slices
        
        warp = cv.CreateMat( 2, 3, cv.CV_32FC1 )
        
        warp[0,0] = 1.0
        warp[0,1] = 0.0
        warp[1,0] = 0.0
        warp[1,1] = 1.0
        warp[0,2] = 0.0
        warp[1,2] = 0.0
        
        volumes = []
        for iz in range(self.nz):
            slice = mrcVolume.getZSliceAt(iz)
            slice = slice[:,::-1].T.copy()
            slice = cv.fromarray(slice)
            slice_ = cv.CreateMat( self.ny, self.nx, cv.CV_32FC1 )
            cv.Convert( slice, slice_ )
            volumes.append(slice_)
        
        
        slice = cv.CreateMat( self.nx, self.ny, cv.CV_32FC1)
        final_slice = cv.CreateMat( self.ny, self.nx, cv.CV_32FC1)
        
        for index,angle in enumerate(self.angles2process):
            
            sys.stdout.write( "Angle: %i: Z-> " %angle )
         
            angle_in_rad = numpy.radians(angle)
        
            src_file = "%s/%s_%.1f_f.mrc" %( self.filt_dir,self.basename,angle )
            fsrc = mrc.MRCFile(src_file)
            
            u = fsrc.getZSliceAt(0)
            u = u[:,::-1].T.copy()
            src_ = cv.fromarray(u)
            
            src = cv.CreateMat( self.ny, self.nx, cv.CV_32FC1 )
            cv.Convert( src_, src )
        
            for iz in range( self.nz ):
                
                Z = self.Zmin + iz
                
                mat = self.stretch_image( src, angle_in_rad, Z=Z, extra_translation=(self.t[index,0],self.t[index,1]) )
                
                sys.stdout.write( '%i..' %Z )
                sys.stdout.flush()
                
                cv.Add( mat, volumes[iz], final_slice )
                cv.Copy( final_slice, volumes[iz] )
                
            sys.stdout.write("\n")
                            
                            
        for iz in range( self.nz ):
            u = numpy.asarray(volumes[iz])
            u = u[::-1,:].T
            mrcVolume.setZSliceAt( iz, u )
            
        mrcVolume.updateHeader()
   
        
    def save(self):
        '''Save the data configuration'''
        
        if len(self.angles2process)==0: return
        
        # All the angles have been processed
        for angle,t in zip(self.angles2process,self.t.tolist()):
            if not angle in self.processed:
                self.processed.append(angle)
                self.processed_t.append(t)

        f = open(self.config_file,'wb')
        f.write("%s: [%.6f,%.6f,%.6f]" %(rotation_axis,self.rot_axis[0],self.rot_axis[1],self.rot_axis[2]))
        
        cPickle.dump((self.processed,self.processed_t),open(self.status_file,'wb'))
        
        if self.showVolume:
            mrcVolume = os.path.join( self.vol_dir, '%s_vol.mrc' %( self.basename ) )
            os.system( 'imod %s' %mrcVolume )
        
        if self.makeAnimatedGif:
            txbr.utilities.make_animated_gif('%s_vol.mrc' %( self.basename ), self.vol_dir)
        
        

if __name__ == '__main__':
    
    directory = "."
    workDirectory = directory
    basename = "HPF3Viewprep070610b"
    extension = TIF
    start_angle = -60.0
    inc_angle = 2.0
    
    #directory = "."
    #basename = "zap"
    #extension = JPG
    #start_angle = -60.0
    #inc_angle = 1.0
    
    #directory = "."
    #basename = "npoA_5d2_tomo"
    #extension = TIF
    #start_angle = -60.0
    #inc_angle = 2.0    
    
    resetNames = False
    reconstructOnTheFly = True
    
    if resetNames:
        reset( basename, extension=extension, start_angle=start_angle, inc_angle=inc_angle )
    
    if reconstructOnTheFly:
        rec = ReconstructionOnTheFly( directory, workDirectory, basename, extension=extension )
        rec.update()
    
