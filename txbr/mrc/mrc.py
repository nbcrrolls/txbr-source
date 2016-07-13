import os
import numpy
import mrcfile
import util

DEFAULT_WIDTH = 100
DEFAULT_HEIGHT = 100
DEFAULT_DEPTH = 1


class MRCFile:
    '''A class to handle basic operations on MRC files'''

    def __init__( self, filename, template=None):

        self.filename = filename
        
        # Init default values
        
        self.nx, self.ny, self.nz = DEFAULT_WIDTH, DEFAULT_HEIGHT, DEFAULT_DEPTH
        self.x_origin, self.y_origin, self.z_origin = 0.0, 0.0, 0.0
        self.description = []
        
        # Load Values from file if possible

        if os.path.exists(filename) and os.path.isfile(filename):
            self.getHeader()

        # Create a header from the template if it exists

        if template!=None and os.path.exists(template):
            mrcfile.create_header_from(filename,template)


    def getHeader( self ):
        '''Read the header informations from the MRC file'''

        self.mode = mrcfile.mode(self.filename)
        self.nx, self.ny, self.nz = mrcfile.size(self.filename)
        self.sx, self.sy, self.sz = self.getScales()
        self.x_origin, self.y_origin, self.z_origin, self.xlen, self.ylen, self.zlen = mrcfile.coords(self.filename)
        self.nxstart, self.nystart, self.nzstart, self.mx, self.my, self.mz = mrcfile.grid(self.filename)
        self.description = mrcfile.description(self.filename)


    def setHeader( self, nx, ny, nz, mode=2, x_origin=0.0, y_origin=0.0, z_origin=0.0, stats=None ):
        '''Set the MRC header file.
        Variables nx, ny, nz: Effective dimensions of the MRC file.
        Variables x_origin, y_origin, z_origin: Real origin Coordinates of the MRC file.
        Variable mode: mode of the MRC file. The default value (2) corresponds to float.
        '''

        self.nx = int(nx)
        self.ny = int(ny)
        self.nz = int(nz)
        
        self.x_origin = x_origin
        self.y_origin = y_origin
        self.z_origin = z_origin
        
        self.mode = int(mode)

        mrcfile.write_header( self.filename, nx, ny, nz, x_origin, y_origin, z_origin, mode )
        
        if stats!=None: 
            mrcfile.update_header(self.filename,stats)


    def updateHeader( self, stats=None ):
        '''stats: a dictionary that could specify m0,mean,m2,min,max'''

        (self.min,self.mean,self.max) = mrcfile.update_header(self.filename,stats)


    def update(self):
        '''Update the MRCFile instance header info from values in the file'''
        
        self.getHeader()
        

    def stats(self):
        '''Get the statistics'''
        
        return mrcfile.stats( self.filename )
        
        
    def getScales(self):
        '''Get the scale (in Angstrom per pixel) for the trhee directions of this MRC file'''
        
        return mrcfile.scale(self.filename)
    
        
    def updateOrigin( self, x_origin=None, y_origin=None, z_origin=None ):
        '''Update the origin'''
        
        if x_origin!=None: self.x_origin = float(x_origin)
        if y_origin!=None: self.y_origin = float(y_origin)
        if z_origin!=None: self.z_origin = float(z_origin)
        
        mrcfile.update_origin( self.filename, self.x_origin, self.y_origin, self.z_origin )
        
        
    def updateLength( self, xlen=None, ylen=None, zlen=None ):
        '''Update the origin'''
        
        if xlen!=None: self.xlen = float(xlen)
        if ylen!=None: self.ylen = float(ylen)
        if zlen!=None: self.zlen = float(zlen)
        
        mrcfile.update_length( self.filename, self.xlen, self.ylen, self.zlen )
        
 
    def updateGrid( self, nxstart=None, nystart=None, nzstart=None, mx=None, my=None, mz=None ):
        '''Update the Grid parameters'''
        
        if nxstart!=None: self.nxstart = int(nxstart)
        if nystart!=None: self.nystart = int(nystart)
        if nzstart!=None: self.nzstart = int(nzstart)
        if mx!=None: self.mx = int(mx)
        if my!=None: self.my = int(my)
        if mz!=None: self.mz = int(mz)
        
        mrcfile.update_grid( self.filename, self.nxstart, self.nystart, self.nzstart, self.mx, self.my, self.mz )

    
    def getXSliceAt( self, slice ):
        '''Returns an numpy array for the X slice at index slice.
        Variable slice: The index of the X slice.'''

        return mrcfile.read_slice( self.filename, 'X', slice )


    def setXSliceAt( self, slice, u ):
        '''Set the X slice at index slice for this MRC file.
        Variable slice: The index of the X slice.
        Variable u: A numpy array that contains the density values for this X slice.
        '''

        mrcfile.write_slice( self.filename, 'X', slice, u)


    def getMiddleSliceInX( self ):
        '''Returns the middle slice in the X direction.
        '''

        return mrcfile.read_slice( self.filename, 'X', self.nx/2 )



    def getYSliceAt( self, slice ):
        '''Returns an numpy array for the Y slice at index slice.
        Variable slice: The index of the Y slice.'''

        return mrcfile.read_slice( self.filename, 'Y', slice )


    def setYSliceAt( self, slice, u ):
        '''Set the Y slice at index slice for this MRC file.
        Variable slice: The index of the Y slice.
        Variable u: A numpy array that contains the density values for this Y slice.
        '''

        mrcfile.write_slice(self.filename,'Y',slice,u)


    def getMiddleSliceInY( self ):
        '''Returns the middle slice in the Y direction.
        '''

        return mrcfile.read_slice( self.filename, 'Y', self.ny/2 )


    def getZSliceAt( self, slice ):
        '''Returns an numpy array for the Z slice at index slice.
        Variable slice: The index of the Z slice.'''

        return mrcfile.read_slice( self.filename, 'Z', slice )


    def setZSliceAt( self, slice, u ):
        '''Set the Z slice at index slice for this MRC file.
        Variable slice: The index of the Z slice.
        Variable u: A numpy array that contains the density values for this Z slice.
        '''
        
        mrcfile.write_slice(self.filename,'Z',slice,u)


    def getMiddleSliceInZ( self ):
        '''Returns the middle slice in the Z direction.
        '''

        return mrcfile.read_slice( self.filename, 'Z', int(self.nz/2.0) )
        
    
    def getDescription( self ):
        '''Returns the description labels of this MRC file.'''
        
        self.update()
        return self.description
    
    
    def setDescription( self, labels ):
        '''Sets the description labels for this MRC file.'''

        mrcfile.update_description(self.filename,labels)
        self.update()
        
    
    def copyDescriptionFrom( self, filename ):
        '''Copy the description labels from a MRC file called filename.'''
        
        mrcfile.copy_description( self.filename, filename )
        self.update()
        

    def addLabel( self, label ):
        '''Add a new label for this MRC file description.'''
        
        mrcfile.add_label( self.filename, label )
        self.update()


    def warp( self, transformations, output=None):
        '''Create a MRC stack whose slices are warped from original data.
        transformations: a list containing the transformation (affine 2D) for each
        slice'''

        if type(transformations).__name__=='dict':
            trfs = []
            for index in range(self.nz):
                if transformations.has_key(index):
                    trfs.append(transformations[index])
                else:
                    trfs.append(None)

        if output==None:
            out = self
        else:
            out = MRCFile(output, template=self.filename)

        w = numpy.zeros((self.nx,self.ny))

        for index,trf2D in enumerate(trfs):
            print "Transformation for tilt %i: %s" %(index,trf2D)
            u = self.getZSliceAt(index)
            #u = 1000*(u - numpy.mean(u))/numpy.std(u)
            if trf2D!=None: u = util.warp( u, trf2D )
            out.setZSliceAt(index,u)

        out.updateHeader()


    def fileinfo(self):
        
        mrcfile.header(self.filename)
        

    def info( self ):
        '''Returns informations from the header of this MRC file'''

        self.update()

        info = 'MRC file Name: %s (mode:%i)' %(self.filename,self.mode)
        info += '\nVolume: (nx,ny,nz)=(%i,%i,%i)' %(self.nx,self.ny,self.nz)
        info += '\nSpecimen Dimensions: Origin=(%.2f,%.2f,%.2f)  Size=(%.2f,%.2f,%.2f)' \
                    %( self.x_origin, self.y_origin, self.z_origin, self.xlen, self.ylen, self.zlen )
        info += '\nGrid to load: Origin=(%i,%i,%i)  Size=(%i,%i,%i)' \
                    %( self.nxstart, self.nystart, self.nzstart, self.mx, self.my, self.mz )
        info += '\nScale: (sx,sy,sz)=(%.1f,%.1f,%.1f)' %(self.sx,self.sy,self.sz)
        if len(self.description)!=0:
            info += '\nDescription:\n%s' %('\n'.join(self.description))
        
        return info
    
    
    def __repr__( self ):
        
        return self.info()
    