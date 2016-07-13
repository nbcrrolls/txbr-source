import sys
import os.path
import math
import time
import numpy
import numpy.random
import numpy.linalg
import scipy.optimize
import pylab
import modl
import txbr
import txbr.align
import txbr.filter
import txbr.setup
import modl.structure
import util

from residualMB import ResidualMB
from residualMBOrtho import ResidualMBOrthogonal
from residualCB import ResidualCB
from residualCBOrtho import ResidualCBOrthogonal

from txbr import log

MAXITER = txbr.MAXIMUM_NUMBER_ITERATION
TOL = txbr.TOL

REMAP_ORDER = txbr.REMAP_ORDER
FILTER_MAGNIFICATION = txbr.FILTER_MAGNIFICATION

X_PADDING = txbr.X_PADDING
Y_PADDING = txbr.Y_PADDING
Z_PADDING = txbr.Z_PADDING

SAVE_PLOT = txbr.SAVE_PLOT

modl.structure.USE_MAYAVI = txbr.USE_MAYAVI

N1_DEFAULT = 0
N2_DEFAULT = 1
N3_DEFAULT = 0
N4_DEFAULT = 0

T1_FLAG = 1
T2_FLAG = 2
T3_FLAG = 4
PHI1_FLAG = 8
PHI2_FLAG = 16
PHI3_FLAG = 32
MULTIPLE_BEAM_FLAG = 64
ORTHOGONAL_FLAG = 128

numpy.set_printoptions(precision=4)
numpy.set_printoptions(linewidth=numpy.Infinity)

def __encode_model__( orthogonal, multiple, constant_phi, constant_t ):

    model = 0

    if orthogonal: model += ORTHOGONAL_FLAG
    if multiple: model += MULTIPLE_BEAM_FLAG
    if constant_phi[0]: model += PHI1_FLAG
    if constant_phi[1]: model += PHI2_FLAG
    if constant_phi[2]: model += PHI3_FLAG
    if constant_t[0]: model += T1_FLAG
    if constant_t[1]: model += T2_FLAG
    if constant_t[2]: model += T3_FLAG

    return model

def __decode_model__( model ):

    orthogonal = (model & ORTHOGONAL_FLAG) == ORTHOGONAL_FLAG
    multiple = (model & MULTIPLE_BEAM_FLAG) == MULTIPLE_BEAM_FLAG
    phi_constant = [ (model & PHI1_FLAG) == PHI1_FLAG, (model & PHI2_FLAG) == PHI2_FLAG, (model & PHI3_FLAG) == PHI3_FLAG ]
    t_constant = [ (model & T1_FLAG) == T1_FLAG, (model & T2_FLAG) == T2_FLAG, (model & T3_FLAG) == T3_FLAG ]

    phi_constant = numpy.asarray(phi_constant)
    t_constant = numpy.asarray(t_constant)

    model_rep = __model_representation__(orthogonal,multiple,phi_constant,t_constant)

    log.info("Beam Model: {0:s}".format(model_rep))

    return ( orthogonal, multiple, phi_constant, t_constant )


def __model_representation__( orthogonal, multiple, phi_constant, t_constant, orders=None, flatten_order=None, verbose=True ):
    
    (n1,n2,n3,n4) = orders==None and (0,1,0,0) or orders
    
    rep = ''

    if verbose:

        if orthogonal: rep += 'Orthogonal'
        else:  rep += 'General'
        
        rep += ', Order {0:}'.format(n2)
    
        if multiple: rep += ', Multiple'
        else:  rep += ', Constant'
    
        if orthogonal: rep += ', phi-fixed [{0:},{1:},{2:}]'.format( phi_constant[0], phi_constant[1], phi_constant[2] )
    
        rep += ', t-fixed [{0:},{1:},{2:}]'.format( t_constant[0], t_constant[1], t_constant[2] )

        if flatten_order!=None: rep += ', flatten: {0:s}'.format(str(flatten_order))
    
    else:
        
        if orthogonal: rep += 'Ortho.'
        else:  rep += 'Gen.'
        
        rep += ', n2={0:}'.format(n2)
    
        if multiple: rep += ', Mult.'
        else:  rep += ', Const.'
    
        if orthogonal: rep += ', phi-fixed [{0:},{1:},{2:}]'.format( phi_constant[0], phi_constant[1], phi_constant[2] )
    
        rep += ', t-fixed [{0:},{1:},{2:}]]'.format( t_constant[0], t_constant[1], t_constant[2] )

        if flatten_order!=None: rep += ', flatten: {0:s}'.format( str(flatten_order) )

    return rep



def __evaluate_orthogonal_from_projection__( P ):
    '''
    Given a projection map, this routine evaluates the corresponding
    orthogonal parameters from its linear part.
    '''

    Px = P[::2,:]
    Py = P[1::2,:]

    tan_angle_3 = Px[:,2]/Px[:,1]
    tan_angle_2 = Px[:,3]/numpy.sqrt(Px[:,1]**2 + Px[:,2]**2)
    tan_angle_1 = Py[:,3]/(Py[:,2]*Px[:,1]-Py[:,1]*Px[:,2])

    angle_1 = numpy.arctan(tan_angle_1)
    angle_2 = numpy.arctan(tan_angle_2)
    angle_3 = numpy.arctan(tan_angle_3)

    # Calculating the angles from the projection might lead to problems especially
    # when the value of the coefficients are close to 0. Sometimes in the case of
    # angle_3 it is better to calculate it from its arccos rather than its arctan

    cos_angle_3 = Px[:,1]/numpy.cos(angle_2)

    cos_angle_3 = numpy.where(cos_angle_3>1.0,1.0,cos_angle_3)
    cos_angle_3 = numpy.where(cos_angle_3<-1.0,-1.0,cos_angle_3)

    angle_3 = numpy.arccos(cos_angle_3)

    return numpy.row_stack((angle_1,angle_2,angle_3))

def __re_sample_coefficients__( coeffs, decimate, order=3 ):
    """
    :param coeffs: A list containing the parameters regress.
    :param decimate:
    :param order: The order of the regression (default value:3).
    """

    regular_indices = numpy.array(range(len(coeffs)))

    decimated_indices = regular_indices[::decimate]
    decimated_coeffs = coeffs[::decimate]

    pfit = numpy.polyfit(decimated_indices,decimated_coeffs,order)
    reg_fn = numpy.poly1d(pfit)

    regular_coeffs = reg_fn(regular_indices)
    regular_coeffs[::decimate] = decimated_coeffs[:]

    return regular_coeffs

class TxBRContourAlign:
    '''A class to handle bundle adjustment'''


    def __init__( self, project, model, n1=0, n2=1, n3=0, n4=0, axis=None, flatten_order=None, \
                  remap_order=REMAP_ORDER, remap_magnification=FILTER_MAGNIFICATION, shortcut=True, \
                  do3dPlot=False, doPlot=False, prepare=True, skip_ncg=False, \
                  iter=MAXITER, blocksize=None, \
                  init=False, applyRelativeTransformAtInit=True, evaluateFromScratch=True, \
                  with_contour=True, input = None, fix_series_projection_map=[], \
                  enforceTangencyConstraint = txbr.ENFORCE_TANGENCY_CONSTRAINT, \
                  structureCharacteristicLength = txbr.STRUCTURE_CHARACTERISTIC_LENGTH, \
                  mosaic=False, decimate=None ):
        """
        :param decimate: if None or 1, all the tilt for each series are considered. Otherwise, only
        consider the tilts every *decimate*.
        """

        log.info("Initialize Alignment object")

        self.project = project
        self.processing = False
        self.extraPlot = False
        self.doPlot = doPlot or SAVE_PLOT  # plots are always saved in png files if doPlot is True
        self.showPlot3D = do3dPlot
        self.decimate = decimate==None and 1 or int(decimate)
        
        self.align_directory = project.align_directory
        self.contour_directory = project.contour_directory
        self.filter_directory = project.filter_directory

        if not os.path.exists( self.align_directory ): os.makedirs(self.align_directory)
        if not os.path.exists( self.filter_directory ): os.makedirs(self.filter_directory)
        if not os.path.exists( self.contour_directory ): os.makedirs(self.contour_directory)

        self.project.validate()

        self.orthogonal, self.multiple, self.phi_constant, self.t_constant = __decode_model__( model )

        self.skip_ncg = skip_ncg

        self.n1 = with_contour and n1 or 1
        self.n2 = with_contour and n2 or 1
        self.n3 = with_contour and n3 or 1
        self.n4 = with_contour and n4 or 0
                    
        log.info('Order ({0:},{1:},{2:},{3:})'.format(n1,n2,n3,n4))

        if isinstance(flatten_order, (list, tuple)):
            self.flatten_order = flatten_order[0]
            if len(flatten_order)>1: self.flatten_order_eff = flatten_order[1]
            else: self.flatten_order_eff = flatten_order[0]
        elif isinstance(flatten_order, int):
            self.flatten_order = flatten_order
            self.flatten_order_eff = None   # Projection map and flattened map will be composed
        else:
            self.flatten_order = None
            self.flatten_order_eff = None

        self.remap_order = remap_order
        self.remap_magnification = remap_magnification
        self.shortcut = shortcut
        self.iter = iter
        self.init = init
        self.applyRelativeTransformAtInit = applyRelativeTransformAtInit
        self.evaluateFromScratch = evaluateFromScratch
        self.with_contour = with_contour
        self.input = input
        self.enforceTangencyConstraint = enforceTangencyConstraint
        self.structureCharacteristicLength = structureCharacteristicLength
        self.fix_series_projection_map = fix_series_projection_map
        self.mosaic = mosaic
        
        if blocksize!=None:
            self.project.reconstruction.sizeOfBlock = blocksize

        log.info('Full Minimization: %s' %(not self.shortcut))

        self.parameter_file = os.path.join(self.align_directory,project.basename +'.ctl')

        if axis!=None:
            for series in project.series: series.rotAxis = axis

        self.u_t = [ numpy.asarray(series.rotAxis) for series in project.series ]

        for u in self.u_t:
            log.info('Reference Axis: [{0:.2f},{1:.2f},{2:.2f}]'.format(u[0],u[1],u[2]))

        self.initialize()

        self.initializeResidual(self.n1,1,self.n3,self.n4,prepare)  # Treat first the residual with a linear projection map

        self.doInitialEstimates()
        
        
    def model(self):
        
        return __encode_model__( self.orthogonal, self.multiple, self.phi_constant, self.t_constant)
        
    
    def __model_representation__( self , n1=None, n2=None, n3=None, n4=None ):
        
        if n1==None: n1 = self.n1
        if n2==None: n2 = self.n2
        if n3==None: n3 = self.n3
        if n4==None: n4 = self.n4

        orders = ( n1, n2, n3, n4 )

        flatten_order = None
        if self.flatten_order!=None: flatten_order = [self.flatten_order]
        if self.flatten_order_eff!=None: flatten_order.append(self.flatten_order_eff)

        return __model_representation__( self.orthogonal, self.multiple, self.phi_constant, self.t_constant, orders=orders, flatten_order=flatten_order )
        

    def initialize( self ):
        """
        Load relevant data (Number of series, tilt and objects) from the project model file.
        """

        nseries = self.project.numberOfSeries()

        self.seriesRange = []    # Index range of a series
        self.referenceExposures = []    # Reference Exposure index vs series index
        self.tlts = []    # Number of exposures vs series index
        self.seriesBasenames = []    # Number of exposures vs series index
        self.fix_tilt_range = []
        self.numberOfTilts = 0    # Total number of exposures
        self.numberOfPatches = 0    # Total number of objects
        self.models = []
        self.points = []
        self.lines = []
        self.surfaces = []
        
        if nseries==0: return;

        for series in self.project.series:
            
            if not series.enabled: continue
            
            offset = self.numberOfTilts

            trackModel = series.getTrackingModel()
            trackModel.loadFromFile(keepOnlyPointStructures=not self.with_contour)
            trackModel.tiltOffset = offset

            self.seriesBasenames.append(series.basename)

            number_of_tilts_in_series = int(1+(series.numberOfExposures()-1)/self.decimate)
            ref_tilt_in_series = int(1+(series.indexOfReferenceExposure()-1)/self.decimate)

            self.referenceExposures.append(offset+ref_tilt_in_series)
            self.seriesRange.append(range(offset,offset+number_of_tilts_in_series))
            self.tlts.append(number_of_tilts_in_series)
            self.numberOfTilts += number_of_tilts_in_series
            
            self.models.append(trackModel)
            
            for extension in self.fix_series_projection_map:
                if serie.basename.endswith(extension):
                    self.fix_tilt_range += range(offset,offset+number_of_tilts_in_series)

        self.numberOfPatches = self.models[0].numberOfObjects()
        
        for model in self.models:
            if model.numberOfObjects()!=self.numberOfPatches:
                raise ValueError, "All series should have the same number of objects!"

        log.info('Micrographs Dimensions: (nx,ny)=(%i,%i)' %(self.project.nx,self.project.ny))
        log.info('Total Number of Tilts: %i (%s)' %(self.numberOfTilts,self.tlts))
        log.info('Number of Patches: %i' %self.numberOfPatches)

        # Eventually make sure that the models contain the same object list

        model = self.models[0]

        for object in model.objects:
            indexOfObject = object.indexOfObject
            if object.isAPoint():
                self.points.append(indexOfObject)
            elif object.isALine():
                self.lines.append(indexOfObject)
            else:
                self.surfaces.append(indexOfObject)
                    
        self.pure_bead_markers = len(self.lines)==0 and len(self.surfaces)==0

        log.info('Number of Point Object Structures: %i' %len(self.points))
        log.debug('  %s' %self.points)
        log.info('Number of Line Object Structures: %i' %len(self.lines))
        log.debug('  %s' %self.lines)
        log.info('Number of Surface Object Structures: %i' %len(self.surfaces))
        log.debug('  %s' %self.surfaces)


    def initializeResidual( self, n1, n2, n3, n4, prep ):
        """
        Initialize the residual object (to derive symbolically the cost function
        expressions, their derivatives and hessians - see bundle adjustment).
        """

        ntilt = self.numberOfTilts
        npatch = self.numberOfPatches

        if not self.orthogonal and not self.multiple:

            self.resid = ResidualCB( ntilt, npatch, n1, n2, n3, n4, fix_b = True, \
                                    prepare=prep, enforceTangencyConstraint = self.enforceTangencyConstraint, \
                                    structureCharacteristicLength = self.structureCharacteristicLength )

        elif not self.orthogonal and self.multiple:

            self.resid = ResidualMB( ntilt, npatch, n1, n2, n3, n4, fix_b=True, \
                                     prepare=prep, enforceTangencyConstraint= self.enforceTangencyConstraint, \
                                     structureCharacteristicLength = self.structureCharacteristicLength )

        elif self.orthogonal and not self.multiple:

            self.resid = ResidualCBOrthogonal( ntilt, npatch, n1, n2, n3, n4, fix_b=True, \
                                               prepare=prep, enforceTangencyConstraint= self.enforceTangencyConstraint, \
                                               structureCharacteristicLength = self.structureCharacteristicLength )

        elif self.orthogonal and self.multiple:

            self.resid = ResidualMBOrthogonal( ntilt, npatch, n1, n2, n3, n4, fix_b=True, \
                                               prepare=prep, enforceTangencyConstraint= self.enforceTangencyConstraint, \
                                               structureCharacteristicLength = self.structureCharacteristicLength )

        self.resid.t_constant = self.t_constant
        self.resid.phi_constant = self.phi_constant


    def doInitialEstimates( self ):
        """
        Initialize estimates for the residual coefficients.
        """
        
        log.info("Do initial estimates!")
        
        if self.project.numberOfSeries()==0: return

        self.tracks, self.mask, self.glomap = self.initializeTracks()

        self.initializeTraces()

        mode = self.init and 'project' or 'ideal'
        set = self.mosaic and 'global' or 'local'
       
        self.angles, self.P = self.initializeProjectionMap( mode=mode, set=set )    # Ideal Projection
       
        self.resid.P = self.P
        self.resid.angles = self.angles
        #self.resid.angles = numpy.degrees(__evaluate_orthogonal_from_projection__(self.P))
        
        self.resid.tracks = self.tracks
        self.resid.Emask = self.mask
        self.resid.points = self.points
        self.resid.lines = self.lines
        self.resid.surfaces = self.surfaces
        self.resid.fixed_origin = numpy.array([ self.project.nx/2.0, self.project.ny/2.0, 0.0 ])

        self.resid.tlts = self.tlts
        
        if self.doPlot and False:
            self.plot_individual_contour_mappings()

        if self.evaluateFromScratch or not os.path.exists(self.parameter_file):

            log.info("Evaluate new parameters from scratch!")

            self.resid.a = self.resid.initializeACoefficients( fromSingleSeries=False )
            self.resid.b = self.resid.initializeBCoefficients()
            self.resid.b_scaling = self.resid.initializeBScalingCoefficients()
            self.resid.c = self.resid.initializeCCoefficients()
            self.resid.d = self.resid.initializeDCoefficients()

        else:    # Read them from a file

            log.info("Load parameters from %s!", self.parameter_file)

            [ self.resid.a, self.resid.b,  self.resid.b_scaling, self.resid.c, self.resid.d ] = \
                        self.loadParameters(self.parameter_file)

        # Locally the projection map is stored in its general form

        self.a = self.resid.unfold_var1(self.resid.a)
        self.b, self.b_scaling = \
            self.resid.unfold_var2_with_scaling(self.resid.b,self.resid.b_scaling)
        self.c = self.resid.unfold_var3(self.resid.c)
        self.d = self.resid.unfold_var4(self.resid.d)
        


    def initializeProjectionMap( self, mode='ideal', set='local' ):
        """
        Initialize the projection map.
        :param mode: If the variable 'mode' is set to 'ideal': the projection map is calculated from the
                                                  experimental angles.
                                     'project': the projection map is the current one of
                                                  the project.
        :param set: If the variable 'set' is set to 'global': the tilt series are part of one global tilt series and
                                                  thus the rotation axis is close to the middle of the
                                                  the whole reconstruction (image shift mosaic)
                                         'local': the tilt series are independent
        """

        log.info("Initialize the projection map!")

        n1 = self.numberOfTilts
        n2 = self.resid.nn2

        # Build the list of angles

        angles_0 = []

        for iseries,series in enumerate(self.project.series):
            if not series.enabled: continue
            angles_0.append(numpy.asarray(series.tiltAngles[::self.decimate]))

        angles_0 = - numpy.concatenate(angles_0)

        # Tilt angle are defined arount a rotation axis (ey usually)

        angles_0 = angles_0[numpy.newaxis,:]
        angles_0 = numpy.insert(angles_0,[0,1],0,axis=0)

        angles = angles_0

        # Build the projection map

        P = numpy.zeros((2*n1,n2))

        if mode=='ideal':
                      
            for iseries,series in enumerate(self.project.series):

                if not series.enabled: continue

                Origin_src = numpy.mean(series.getSourceFrame(full=False),axis=1)
                Origin_src = numpy.append(Origin_src,0.0)

                if set=='local':
                    Origin_dest = Origin_src.copy()
                else:
                    Origin_dest = numpy.mean(self.project.reconstruction.getFrame(full=False),axis=1)
                    Origin_src = Origin_src + Origin_dest - numpy.mean(series.getDestinationFrame(full=False),axis=1)
                
                u = self.u_t[iseries]    # Rotation axis


#                print "Frame"
#                print self.project.reconstruction.getFrame()
#                print "Origin: (src)->%s   (dest)->%s" %(str(Origin_src),str(Origin_dest))
#                print "u=%s   %f" %(str(u),numpy.degrees(numpy.arctan(u[1]/u[0])))
                
                x_range = 2*numpy.array(self.seriesRange[iseries])
                y_range = 2*numpy.array(self.seriesRange[iseries]) + 1

                cos_phi = numpy.cos(angles[1,self.seriesRange[iseries]]*math.pi/180.0)
                sin_phi = numpy.sin(angles[1,self.seriesRange[iseries]]*math.pi/180.0)

                P[x_range,1] = u[0]*u[0]*(1-cos_phi) + cos_phi
                P[x_range,2] = u[0]*u[1]*(1-cos_phi) - u[2]*sin_phi
                P[x_range,3] = u[0]*u[2]*(1-cos_phi) + u[1]*sin_phi

                P[y_range,1] = u[1]*u[0]*(1-cos_phi) + u[2]*sin_phi
                P[y_range,2] = u[1]*u[1]*(1-cos_phi) + cos_phi
                P[y_range,3] = u[1]*u[2]*(1-cos_phi) - u[0]*sin_phi

                P[x_range,0] = Origin_dest[0] - numpy.dot(P[x_range,1:4],Origin_src)
                P[y_range,0] = Origin_dest[1] - numpy.dot(P[y_range,1:4],Origin_src)
                

                # Handle an eventual transfer rotation

#                (sample_tf_t0,sample_tf_m) = series.getSampleTransferCoordinatesElements()
#                
#                log.info('Reference Transformation for series #%i' %iseries)
#                log.info('   %s' %str(numpy.asarray(sample_tf_t0)))
#                log.info('   %s' %str(numpy.asarray(sample_tf_m)))
#
#                sample_tf_m = numpy.asarray(sample_tf_m).T.tolist()
#
#                P[x_range,0] += numpy.dot(P[x_range,1:4],Origin) - numpy.dot(P[x_range,1:4],numpy.dot(sample_tf_m,Origin))
#                P[y_range,0] += numpy.dot(P[y_range,1:4],Origin) - numpy.dot(P[y_range,1:4],numpy.dot(sample_tf_m,Origin))
#
#                P[x_range,1:4] = numpy.dot(P[x_range,1:4],sample_tf_m)
#                P[y_range,1:4] = numpy.dot(P[y_range,1:4],sample_tf_m)
                
#                relativeTransform = series.getRelativeReferenceTransform()
#                
#                P[x_range,0] += numpy.dot(P[x_range,1:4],relativeTransform.T)
#                P[y_range,0] += numpy.dot(P[y_range,1:4],relativeTransform.T)
#
#                P[x_range,1:4] = numpy.dot(P[x_range,1:4],relativeTransform.M)
#                P[y_range,1:4] = numpy.dot(P[y_range,1:4],relativeTransform.M)                
            
        elif mode=='project':
            
            # Make the distinction if the projection map were evalated as a whole or
            # individually

            offset = 0

            for series in self.project.series:

                if not series.enabled: continue
                
                nproj = min(n2,series.projection.numberOfTerms)
                nexposure = int(1+(series.projection.numberOfExposures()-1)/self.decimate)

                for iexp in range(nexposure):
                    for iproj in range(nproj):
                        tilt_index = iexp*self.decimate
                        P[offset+2*iexp,iproj] = series.projection.x_coefficients[tilt_index][iproj]
                        P[offset+2*iexp+1,iproj] = series.projection.y_coefficients[tilt_index][iproj]

                offset +=2*nexposure
                
        # Apply relative transform between the series. The projection should already be decimated
        # in the tilt angle direction
        
        if self.applyRelativeTransformAtInit:
            
            for iseries,series in enumerate(self.project.series):
                
                if not series.enabled: continue
                
                log.info("Apply relative transform at initialization for series %s" %series.basename)
            
                x_range = 2*numpy.array(self.seriesRange[iseries])
                y_range = 2*numpy.array(self.seriesRange[iseries]) + 1

                relativeTransform = series.getRelativeReferenceTransform()
                
                P[x_range,0] += numpy.dot(P[x_range,1:4],relativeTransform.T)
                P[y_range,0] += numpy.dot(P[y_range,1:4],relativeTransform.T)

                P[x_range,1:4] = numpy.dot(P[x_range,1:4],relativeTransform.M)
                P[y_range,1:4] = numpy.dot(P[y_range,1:4],relativeTransform.M)             
                

        # Recalculate the angles from the projection...
        # How is it different from angles_0

        angles = numpy.degrees(__evaluate_orthogonal_from_projection__( P ))

        return angles, P

        

    def initializeTracks( self ):
        """
        This routine initializes the positions of the marker objects.
        The track object contains the individual polynomial approximation for the contour traces.
        """

        # Establish first for which tilt index a tracking object is defined.
        # This is the mask

        log.debug('Calculate the mask. ')

        mask = numpy.zeros((self.numberOfTilts,self.numberOfPatches))

        for model in self.models:
            for object in model.objects:
                iobj = object.indexOfObject
                indicesOfViews = [ int(model.tiltOffset+index/self.decimate) for index in object.indicesOfViews() if index%self.decimate==0 ]
                mask[indicesOfViews,iobj] = 1
                
        # Generate polynomial approximations on the objects
        # This will be used as an initial estiates for building the structure
        # surfaces

        log.debug('Calculate the globlal polynomial approximation for the contours. ')

        n3 = self.resid.n3
        nn3 = self.resid.nn3

        glomap = []

        for iseries,model in enumerate(self.models):
            objectMapping = model.polynomialMapping(n3,u_t=self.u_t[iseries],doPlot=self.extraPlot)
            for mapping in objectMapping:
                px,py,indicesOfViews = mapping
                indicesOfViews = [ model.tiltOffset+index/self.decimate for index in indicesOfViews if index%self.decimate==0 ]
                skipViews = numpy.array([True for i in range(self.numberOfTilts)])
                skipViews[indicesOfViews] = False
                glomap.append((px,py,skipViews,self.u_t,'2D'))
                #glomap.append((px,py,skipViews,self.u_t,'3D'))

        # Generates individual mapping for each contours

        log.debug('Generate individual contour mappings and create the tracks object.')

        tracks = numpy.zeros((self.numberOfTilts,self.numberOfPatches,2,nn3))

        for iseries,model in enumerate(self.models):
            for object in model.objects:
                iobject = object.indexOfObject
                for contour in object.contours:
                    (pc_x,pc_y) = contour.polynomialMapping(n3,u_t=self.u_t[iseries],doPlot=self.extraPlot)
                    z = int(round(contour.zExtensionSet().pop(),0))
                    if int(z%self.decimate)!=0: continue
                    itilt = int(model.tiltOffset + z/self.decimate)
                    tracks[itilt,iobject,0,:nn3] = pc_x==None and numpy.nan or pc_x.coeff[:]
                    tracks[itilt,iobject,1,:nn3] = pc_y==None and numpy.nan or pc_y.coeff[:]

        return ( tracks, mask, glomap )


    def initializeTraces(self):
        """
        Calculate the trace angles and deduce the orientation of reference axes.
        This operation is meaningful, when the tilt series is roughly aligned (from
        cross-correlation routines). If it is not the case, the rotation axis direction
        won't be necessarly accurate.
        """

        nseries = self.project.numberOfEnabledSeries()

        log.info('Calculate Reference Rotation Axis from Trace Angles...')

        if self.tracks.size==0: return

        pointTracks = self.tracks[:,self.points,:,:1]
        pointMask = self.mask[:,self.points]

        pointTracks = pointTracks.reshape(pointTracks.shape[:3])

        for iseries in range(nseries):

            rg = self.seriesRange[iseries]
            (ref_angle,angle) = txbr.align.trace_angle(pointTracks[rg],pointMask[rg],self.u_t[iseries])
            
            self.project.getSeries(self.seriesBasenames[iseries]).setResidualTraceAngle(angle)
            
            angle_0 = self.project.getSeries(self.seriesBasenames[0]).residualTraceAngle
            angle_s = self.project.getSeries(self.seriesBasenames[iseries]).residualTraceAngle
            dev = angle_s - angle_0
            if (abs(dev)>numpy.pi/2.0) :
                angle = numpy.pi + angle_s
                self.project.getSeries(self.seriesBasenames[iseries]).setResidualTraceAngle(angle)
            
            log.info('Reference Angle: %f rad / %f deg (Series #%i)' %(ref_angle,numpy.degrees(ref_angle),iseries))
            log.info('Residual Track Angle: %f rad / %f deg (Series #%i)' %(angle,numpy.degrees(angle),iseries))

            log.info('Update the rotation axis...')

            log.info('Old rotation axis: %s' %str(self.u_t[iseries]))

            c = numpy.cos(-angle)
            s = numpy.sin(-angle)
#            c = numpy.cos(angle)
#            s = numpy.sin(angle)
            R = numpy.array([[c,s,0.0],[-s,c,0.0],[0.0,0.0,1.0]])


            self.u_t[iseries] = numpy.dot(R,self.u_t[iseries])
            #self.u_t = numpy.dot(R,[0.0,1.0,0.0])

            self.project.series[iseries].rotAxis = self.u_t[iseries]
            self.u_t[iseries] = self.project.series[iseries].rotAxis

            # TEST
            self.project.reconstruction.rotAxis = self.u_t[iseries]
            # TEST

            #print self.project.series[iseries].rotAxis

            filtering_angle = self.project.series[iseries].getFilteringAngle()

            log.info('New rotation axis:  %s' %str(self.u_t[iseries]))
            log.info('Filtering Angle:    %f rad    %f deg' %(filtering_angle,numpy.degrees(filtering_angle)))

        # ############## 'TO UPDATE'  #####################


        log.info('Calculate relative rotation angle between series...')
        
        if pointTracks.size!=0:

            try:

                relativeReferenceTransforms, relativeOrientations = \
                                    txbr.align.relative_rot_transform( pointTracks[self.referenceExposures], \
                                                                       pointMask[self.referenceExposures], \
                                                                       index_ref_series=0 )      
                                                    
                for index in range(len(relativeReferenceTransforms)):
                    if (relativeReferenceTransforms[index]!=None): continue
                    relativeReferenceTransforms[index] = util.AffineTransformation3D(numpy.eye(3,4,1))
                          
                # Check the sign
                relativeOrientations = [ -r.angles() for r in relativeReferenceTransforms ]

                for iseries,seriesName in enumerate(self.seriesBasenames):
                    
                    series = self.project.getSeries(seriesName)
                    
                    series.setSampleOrientation(relativeOrientations[iseries])
                    series.setRelativeReferenceTransform(relativeReferenceTransforms[iseries])
                    
                    orientations = relativeOrientations[iseries]

                    log.info( 'Relative Angle: %f rad / %f deg (Series #%i/Series #%i)' \
                              %(orientations[2],numpy.degrees(orientations[2]),iseries,0) )
                    
                    log.info( 'Relative Transformation (Series #%i/Series #%i):' %(iseries,0) )
                    
                    log.info( '%s' %relativeReferenceTransforms[iseries] )
                
            except:
                
                log.warning("Unexpected error:", sys.exc_info()[0])
                log.warning('Skipped calculation of relative rotation angle between series...')
                
                pass


    def initializeSkipping( self, _fix_a, _fix_b, _fix_c, _fix_d, projection=None, fix_tilt_range=[] ):
        '''Determine which variables should be kept constant'''

        self.resid.skip_a = self.resid.initializeASkipping()
        self.resid.skip_b = self.resid.initializeBSkipping(projection,fix_tilt_range)
        self.resid.skip_c = self.resid.initializeCSkipping()
        self.resid.skip_d = self.resid.initializeDSkipping()

        skip = []

        if not _fix_a: skip += self.resid.skip_a
        if not _fix_b: skip += self.resid.skip_b
        if not _fix_c: skip += self.resid.skip_c
        if not _fix_d: skip += self.resid.skip_d

        return numpy.asarray(skip)


    def extract3DPointMarkers(self):
        '''Extract the XYZ locations of the point-like structures. We also make sure
        that points that do not have any mark on any tilt are removed as well.'''
        
        XYZ = numpy.resize(self.a,(3,self.resid.nn1,self.numberOfPatches))
        
        MINIMUM_TILT_NUMBER = 5
        
#        XYZ = XYZ[:,0,self.points]
#
#        mask = numpy.sum(self.mask[:,self.points],axis=0)
#        mask = numpy.where(mask>MINIMUM_TILT_NUMBER)[0]
#        #mask = numpy.any(self.mask[:,self.points],axis=0)


        obj_mask = numpy.zeros(self.numberOfPatches)
        obj_mask[self.points] = 1
        obj_mask = obj_mask*numpy.sum(self.mask[:,:],axis=0)
        obj_mask = numpy.where(obj_mask>MINIMUM_TILT_NUMBER)[0]
        
        XYZ = XYZ[:,0,obj_mask]
        
        return XYZ.transpose(),obj_mask
    

    def process( self, structureToInitialize=None):
        '''Process the alignment.'''
        
        if self.project.numberOfSeries()==0: 
            log.warning('No Series to process!')
            return

        if self.tracks.size==0: # fiducialless case
            self.reFrameSeries()
            self.saveProjectionMapToProject()
            return
        
        # Eventually reset the source input of the project
        
        if self.input!=None:
            self.project.reconstruction.input = self.input
            for series in self.project.series:
                series.input = self.input
            
        if structureToInitialize=="-":
            structureToInitialize = range(self.numberOfPatches)
        elif structureToInitialize!=None:
            structureToInitialize.extend(self.points)
        
        self.initializeStructureSet(structureToInitialize)
        
        if self.doPlot and False:
            self.plot_contour_reprojection()
            self.plot_Initial_Estimates_Reprojection()
            
        if structureToInitialize!=None: 
            return

        log.info('Start the optimization!')

        _fix_a = False
        _fix_b = False
        _fix_c = True
        _fix_d = False
        
        fix_tilt_range = self.fix_tilt_range

        if len(self.lines)==0 and len(self.surfaces)==0:
            _fix_d = True    # To be faster

        self.resid.fixVariables(fix_a=_fix_a,fix_b=_fix_b,fix_c=_fix_c,fix_d=_fix_d)

        skip = self.initializeSkipping(_fix_a,_fix_b,_fix_c,_fix_d,projection='linear',fix_tilt_range=fix_tilt_range)

        var = self.resid.joinParameters(_fix_a,_fix_b,_fix_c,_fix_d) # We keep a copy of the initial variables

        self.processing = True
                
        log.info('Bundle Adjustment with Contour Alignment Processing...')

        xopt = []    # Variables used for the minimization
                    
        if not _fix_a:
            self.resid.a = self.reCenterTracks()
            xopt += self.resid.a
        if not _fix_b:
            # Relax the projection map if the beam model is not orthogonal
            if not isinstance(self.resid,ResidualMBOrthogonal) and not isinstance(self.resid,ResidualCBOrthogonal):
                self.resid.b = self.resid.relaxProjection(order=(0,1))
            else:
                b = self.resid.relaxProjection(order=(0,1))
            
                #b = self.resid.unfold_var2(b)
                print 'len(b)=%i    2*self.numberOfTilts*self.resid.nn2=%i    %i' %(len(b),2*self.numberOfTilts*self.resid.nn2,2*self.numberOfTilts*4)
                P = numpy.array(b)
                
                P.resize((2,self.resid.nn2,self.numberOfTilts))
                P = numpy.rollaxis(P,2,0)
                P = P[:,:,:4]
                P = numpy.resize(P,(2*self.numberOfTilts,4))
                
                angl = __evaluate_orthogonal_from_projection__(self.P)
                
                angl = numpy.ravel(angl)
                angl_ = angl.tolist()
                
                self.resid.b = angl_ + self.resid.b[3*self.numberOfTilts:]

            xopt += self.resid.b
            
        if not _fix_c: xopt += self.resid.c
        if not _fix_d: xopt += self.resid.d
                
        xopt = numpy.array(xopt)

        if not self.skip_ncg:
            xopt = scipy.optimize.fmin_ncg( self.resid.error_struct, xopt, self.resid.der_error, args=(skip,),
                                        fhess=self.resid.hess_error, avextol=TOL, maxiter=self.iter )

        log.info('xopt=' + str(xopt)[:60] + ' ...')

        xopt = numpy.where(skip,var,xopt)

        self.resid.storeParameters(xopt,_fix_a,_fix_b,_fix_c,_fix_d)

        # Save order-1 alignment after a quick sample reorientation

        self.a = self.resid.unfold_var1(self.resid.a)
        self.b = self.resid.unfold_var2(self.resid.b)
        self.b, self.b_scaling = \
                self.resid.unfold_var2_with_scaling(self.resid.b,self.resid.b_scaling)
        self.c = self.resid.unfold_var3(self.resid.c)
        self.d = self.resid.unfold_var4(self.resid.d)

        self.reorient()

        self.saveProjectionMapToProject( n2=1, extension='.order-1')

        # Take care of the remaining non linear coefficients

        if self.n2>self.resid.n2 and not self.shortcut:

            log.info('Processing high order coefficients...')

            self.resid.extendProjectionMapOrder(self.n2)

            var = self.resid.joinParameters(_fix_a,_fix_b,_fix_c,_fix_d)

            skip = self.initializeSkipping(_fix_a,_fix_b,_fix_c,_fix_d,projection='non-linear',fix_tilt_range=fix_tilt_range)

            xopt = var.copy()

            xopt = scipy.optimize.fmin_ncg( self.resid.error_struct, xopt, self.resid.der_error, args=(skip,),
                                        fhess=self.resid.hess_error, avextol=TOL, maxiter=self.iter )

            log.info('xopt=' + str(xopt)[:60] + ' ...')

            xopt = numpy.where(skip,var,xopt)

            self.resid.storeParameters( xopt, _fix_a, _fix_b, _fix_c, _fix_d )

           # self.resid.b = self.resid.relaxProjection()
           
        if self.n2>self.resid.n2 and self.shortcut:
            
            if self.multiple:
                full_reset = True
            else:
                full_reset = False

            self.resid.extendProjectionMapOrder(self.n2, full_reset=full_reset)

            skip = self.initializeSkipping(_fix_a,_fix_b,_fix_c,_fix_d,projection='non-linear',fix_tilt_range=fix_tilt_range)

            self.resid.b = self.resid.relaxProjection()
            
            (self.resid.b,self.resid.b_scaling) = self.resid.evaluateScalingParameters()

        # Unfold the data into the general reference frame

        self.a = self.resid.unfold_var1(self.resid.a)
        self.b = self.resid.unfold_var2(self.resid.b)
        self.b, self.b_scaling = \
                self.resid.unfold_var2_with_scaling(self.resid.b,self.resid.b_scaling)
        self.c = self.resid.unfold_var3(self.resid.c)
        self.d = self.resid.unfold_var4(self.resid.d)

        log.info('len(a)=%i  len(b)=%i  len(c)=%i  len(d)=%i' %(len(self.a),len(self.b),len(self.c),len(self.d)))
        log.info('a=' + str(self.a)[:60] + ' ...')
        log.info('b=' + str(self.b)[:60] + ' ...')
        log.info('c=' + str(self.c)[:60] + ' ...')
        log.info('d=' + str(self.d)[:60] + ' ...')

        if self.flatten_order!=None:

            self.flatten()    # Effective order of the projection maps can be increased

        else:
        
            self.reorient()
            
        self.info()

        self.solveGaugeAmbiguity()
        
        self.evaluateRotationAxisDirectly()
            
        self.info()

        self.generateAlignLogFile()
        
        self.evaluateAlignmentTransform()
        
        self.reFrameSeries()
        
        self.evaluateRemapTransform()
        
        self.saveStructureSet()

        self.saveProjectionMapToProject()
        
        self.saveParameters(self.parameter_file)
        
#       print self.resid.getMagnificationParameters()

        if self.doPlot:
            
            self.plot_projections()
            
            self.plot_angles()
            
            self.plot_magnifications()

   #         self.plot_contour_reprojection()

        if self.doPlot and self.extraPlot:

            self.plot_contour_reprojection()
            
            for ipatch in self.surfaces:
                self.plot_contour_reprojection(ipatch=ipatch)

        self.processing = False


    def reCenterTracks(self):
        '''This routine re-center the calculated 3D track positions in the camera plane
        relatively to their original projected values.'''

        nn1 = self.resid.nn1
        nobj = self.numberOfPatches

        tilt_ref = self.referenceExposures[0]
        mask_ref = self.mask[tilt_ref,:]

        Xmarks = self.tracks[tilt_ref,:,0,0]*mask_ref
        Ymarks = self.tracks[tilt_ref,:,1,0]*mask_ref

        Cix = numpy.mean(Xmarks)
        Ciy = numpy.mean(Ymarks)

        A = numpy.resize(self.a,(3,nn1,nobj))
        X = A[0,0,:]
        Y = A[1,0,:]

        Cfx = numpy.mean(X*mask_ref)
        Cfy = numpy.mean(Y*mask_ref)

        A[0,0,:] = A[0,0,:] - Cfx + Cix
        A[1,0,:] = A[1,0,:] - Cfy + Ciy

        return A.ravel().tolist()


    def evaluateRotationAxisDirectly( self ):

        log.warning('Directly evaluate the rotation axis')

        nterm = self.resid.nn2
        ntilt = self.numberOfTilts

        B = numpy.resize(self.b,(2,nterm,ntilt))

        self.u_t = []

        for index,basename in enumerate(self.seriesBasenames):

            log.debug('Serie %i: %s' %(index,basename))

            range = self.seriesRange[index]

#            for itilt in range:
#                eigenvalues,eigenvectors = numpy.linalg.eig(B[:,1:3,itilt])
#                print "%.2f: %s    %.2f: %s" %(eigenvalues[0],str(eigenvectors[0]),eigenvalues[1],str(eigenvectors[1]))
#                if itilt==0: N = ( eigenvectors[0][0], eigenvectors[0][1], 0.0 )
#
#            N = numpy.array(N)
#
#            log.info('Serie %i: %s -> %s' %(index,basename,str(N)))

            itiltref = self.referenceExposures[index]

            Mref = numpy.row_stack((B[0,1:4,itiltref], B[1,1:4,itiltref],(0.0,0.0,1.0)))
            Mref = scipy.linalg.inv(Mref)

            M = numpy.zeros((2*len(range),3))


#            print ntilt
#            print B.shape

            for index,itilt in enumerate(range):
#                print itilt
                M[2*index,:] = numpy.tensordot(Mref,B[0,1:4,itilt],axes=((0),(0)))
                M[2*index,0] -= 1
                M[2*index+1,:] = numpy.tensordot(Mref,B[1,1:4,itilt],axes=((0),(0)))
                M[2*index+1,1] -= 1

            U,s,Vh = scipy.linalg.svd(M)
            N = Vh[2,:]/numpy.linalg.norm(Vh[2,:])

            log.warning("Evaluate the rotation axis of series #%i (%s) to : %s" %(index,basename,N))


            serie = self.project.getSeries(basename)
            serie.rotAxis = N

            self.u_t.append(N)
            self.project.reconstruction.rotAxis = N



    def solveGaugeAmbiguity( self ):
        '''Gauge Ambiguity Routine. The projection maps should be corrected to become
        the best orthogonal solution. This should be done using geometrical consideration.'''

        self.reFrameSeries()
        
        recFrame = self.project.reconstruction.getFrame(full=False)

        log.info('Solve Gauge ambiguity')

        nterm = self.resid.nn2
        ntilt = self.numberOfTilts

        B = numpy.resize(self.b,(2,nterm,ntilt))
        
        self.M0 = []
        self.u_t = []

        for index,basename in enumerate(self.seriesBasenames):

            log.debug('Serie %i: %s' %(index,basename))

            serie = self.project.getSeries(basename)
            range = self.seriesRange[index]
            
            src_frame = serie.getSourceFrame(full=False)
            dest_frame = serie.getDestinationFrame(full=False)
                          
            ref_t_elements = serie.getTransferCoordinatesElements()
            smpl_t_elements = serie.getSampleTransferCoordinatesElements()
            
            ref_transform = util.AffineTransformation3D(numpy.column_stack(ref_t_elements))
            smpl_transform =  serie.getRelativeReferenceTransform()

#            print B[:,:4,range]
#            print src_frame
#            print dest_frame
#            print ref_transform
#            print smpl_traprintnsform
#            print serie.basename
#            print self.model()

            try:

                M0,N = txbr.align.calc_rotation_axis( B[:,:4,range], src_frame, dest_frame, ref_transform, smpl_transform, \
                                                   doPlot=self.extraPlot, directory=self.align_directory, \
                                                   basename=serie.basename, model=self.model() )
            except:

                log.error("Problem in the solving the gauge for series %s" %basename)

                M0 = [0.0,0.0,0.0]
                N = self.project.reconstruction.rotAxis
            
            self.project.reconstruction.rotAxis = N
            
            self.M0.append(M0)
            self.u_t.append(N)
            
            log.info("Series %s: M0=[%.1f,%.1f,%.1f]   u=[%.2f,%.2f,%.2f]" 
                     %(basename, M0[0], M0[1], M0[2], N[0], N[1], N[2]))
        

    def reorient(self):
        '''This routine uses the point-like marker (gold fiducial) positions to reorient 
        the reconstructed volume so its normal axis correspond to ez.
        Both the 3D locations of the markers and the projection maps will be changed adequatly
        (similar process to the gauge ambiguity solution).'''
        
        log.info("Evaluate reorientation for the 3D volume...")
        
        if len(self.points)==0: 
            log.info("No point-like markers to reorient the volume")
            return

        nobj = self.numberOfPatches
        ntilt = self.numberOfTilts
        nn1 = self.resid.nn1
        
        nx,ny = self.project.nx,self.project.ny

        # First thing, rotate everything for the fiducials to be in the XY plane
        
        XYZ,ojb_mask = self.extract3DPointMarkers()
        
#        print XYZ
#        
#        scale = numpy.max(XYZ,axis=0)
#        txbr.utilities.plot3DMarkers( XYZ[:,0],XYZ[:,1],XYZ[:,2], scale )
        
        
        

        log.info('Reorient specimen with %i point markers.' %XYZ.shape[0])

        bottom_plane, top_plane = txbr.setup.eval_bounds_2(XYZ,XY_range=(self.project.nx,self.project.ny))

        bottom_plane = numpy.asarray(bottom_plane)
        top_plane = numpy.asarray(top_plane)
        medium_plane = numpy.average(numpy.row_stack((bottom_plane,top_plane)),axis=0)

        v = numpy.asarray(medium_plane[1:])    # Normal vector to the specimen surface
        v = v/numpy.sqrt(numpy.dot(v,v))
        
        ez = numpy.array([0.0,0.0,1.0])

        log.info('Bottom plane: %.2f + %.2f x + %.2f y + z = 0' %(bottom_plane[0],bottom_plane[1],bottom_plane[2]))
        log.info('   Top plane: %.2f + %.2f x + %.2f y + z = 0' %(top_plane[0],top_plane[1],top_plane[2]))
        log.info('Medium plane: %.2f + %.2f x + %.2f y + z = 0' %(medium_plane[0],medium_plane[1],medium_plane[2]))

        # Evaluate the reorientation rotation 
        
        M = numpy.zeros((3))
        M[0] = nx/2.0
        M[1] = ny/2.0
        M[2] = -(medium_plane[0]+numpy.dot(medium_plane[1:3],M[:2]))/medium_plane[3]

        M_img = numpy.zeros((3))
        M_img[:2] = M[:2]
        
        center = numpy.array([0.0,0.0,0.0]) # dummy value; will be calculated from M-M_img
        axis = numpy.cross(ez,v)
        theta = numpy.arccos(numpy.dot(ez,v))

        rotation = util.Rotation3D( center, axis, -theta)
        rotation.resetCenter(M,M_img)
        
        log.info("Reorientation Transformation: %s" %rotation)

        # Rotate a and the tracks
        
        A = numpy.resize(self.a,(3,nn1*nobj))
        A = A.transpose()
        A = rotation.forward(A)
        A = A.transpose()
        
        XYZr = rotation.forward(XYZ)  # repeat of rotation A but only on the point-like markers
        
        log.info('Minimum values of a: [%.1f,%.1f,%.1f]' %tuple(numpy.min(XYZr,axis=0)))
        log.info('Maximum values of a: [%.1f,%.1f,%.1f]' %tuple(numpy.max(XYZr,axis=0)))

        # Calculate the new range in Z
        
        calculate_z_range_from_markers = True
        
        if calculate_z_range_from_markers:  # look directly at the rotated markers
            
            XYZmin = numpy.min(XYZr,axis=0)
            XYZmax = numpy.max(XYZr,axis=0)
            
            XYZmin = XYZmin.astype('int')
            XYZmax = XYZmax.astype('int')
            
            # for the X and Y direction we take the default boundary
            
            XYZmin[:2] = 1
            XYZmax[:2] = [nx,ny]
            
        else:
            
            h = numpy.ones((4))
            h[1:] = rotation.center[:]
            
            dmin = numpy.dot(bottom_plane,h)
            dmax = numpy.dot(top_plane,h)
            
            log.info('Distance of center from planes: dmin=%.3f   dmax=%.3f' %(dmin,dmax))
    
            XYZmin = rotation.forward(-rotation.center+dmin*v)
            XYZmax = rotation.forward(-rotation.center+dmax*v)
    
            XYZmin = XYZmin.astype('int')
            XYZmax = XYZmax.astype('int')
    
            XYZmin[:2] = 1
            XYZmax[:2] = [nx,ny]
            
            if XYZmax[2]<=XYZmin[2]:
                log.error('Error in the Z boundarues: Zmin=%.3f   Zmax=%.3f' %(Zmin[2],Zmax[2]))
        
        
        XYZmin -= numpy.array([X_PADDING,Y_PADDING,Z_PADDING])
        XYZmax += numpy.array([X_PADDING,Y_PADDING,Z_PADDING])
        
        # Update the reconstruction object
                
        self.project.reconstruction.setOrigin(*XYZmin)
        self.project.reconstruction.setEnd(*XYZmax)
        self.project.reconstruction.setBottomPlane(0.0,0.0,0.0,XYZmin[2])
        self.project.reconstruction.setTopPlane(0.0,0.0,0.0,XYZmax[2])

        log.info('Reconstruction boundaries: %s -> %s' %(str(XYZmin[2]),str(XYZmax[2])))

        # Now compose the polynomial projection with the rotation

        log.info('Rotate projection map adequatly...')
        
        B = numpy.resize(self.b,(2,self.resid.nn2,ntilt))
        
        Rinv = rotation.inv()
        
        for itilt in range(ntilt):
        
            proj_map = util.PolynomialTransformation3D(self.resid.n2,B[:,:,itilt])
            proj_map = proj_map.compose(Rinv)
            
            B[:,:,itilt] = proj_map.coeffs()[:2,:]
        
        
        log.info('Map reorientation done...')

        self.a = A.ravel().tolist()
        self.b = B.ravel().tolist()


    def flatten(self):
        '''Use the point track markers to modify the projection map to
        be the best orthogonal ones.'''

        if self.flatten_order==None:
            return
        
        log.info("Apply the flattening processing...")
        
        if len(self.points)==0: 
            log.info("No point-like markers to flatten the volume")
            return

        nobj = self.numberOfPatches
        ntilt = self.numberOfTilts
        nn1 = self.resid.nn1
        
        XYZ,obj_mask = self.extract3DPointMarkers()

        log.info('Flatten specimen with %i point markers.' %XYZ.shape[0])

        # Evaluate a polynomial function to represent the boundary surfaces

        center = [ self.project.nx/2.0, self.project.ny/2.0, 0.0]
        
        warp_to_flat_tf, flat_to_warp_tf = txbr.setup.evalFlatteningTransformation( XYZ, self.flatten_order, center, doPlot=self.showPlot3D )

        XYZr = warp_to_flat_tf.forward(XYZ)

        # Evaluate the boundaries from the fiducials
        
        XYZ_min = numpy.min(XYZr,axis=0)
        XYZ_max = numpy.max(XYZr,axis=0)
    
        XYZ_min = XYZ_min.astype('int')
        XYZ_max = XYZ_max.astype('int')
        
        # for the X and Y direction we take the default boundary
        
        XYZ_min[:2] = 1
        XYZ_max[:2] = [self.project.nx,self.project.ny]
        
        XYZ_min -= numpy.array([X_PADDING,Y_PADDING,Z_PADDING])
        XYZ_max += numpy.array([X_PADDING,Y_PADDING,Z_PADDING])
        
        XYZ_min = numpy.rint(XYZ_min)
        XYZ_max = numpy.rint(XYZ_max)
        
        # Update the reconstruction object
                
        self.project.reconstruction.setOrigin(*XYZ_min.tolist())
        self.project.reconstruction.setEnd(*XYZ_max.tolist())
        self.project.reconstruction.setBottomPlane(0.0,0.0,0.0,XYZ_min[2])
        self.project.reconstruction.setTopPlane(0.0,0.0,0.0,XYZ_max[2])   
        
        # Recalculate a and the tracks
        # This will not work with contour alignment
                
        A = numpy.zeros(((3,nn1,nobj)))
        A[:,0,obj_mask] = XYZr.T

        self.a = A.ravel().tolist() # reassign a
                
        # Correct the modified Projection Map
        # Two routes are possible
        # (i) composing the projection map with the inverse of the flattening transform
        # (ii) Relaxing the projection map while the gold markers are in their "flattening" position
        # In the case of the mosaic, since the flattening transform is a global function on whole
        # the tiles, it might be better to effectively take a polynomial "tile flattening function"
        # with a lower order than the global flattening order. A compromise between
        # computational faisability and local flattening has to be made.
            
        log.info('Recalculate the projection map adequatly...')

        #doCompose = self.flatten_order_eff==None or self.flatten_order==self.flatten_order_eff
        doCompose = self.flatten_order_eff==None

        if doCompose:

            B = numpy.resize(self.b,(2,self.resid.nn2,ntilt))

            n2_new = self.flatten_order*self.resid.n2
            nn2_new = util.numberOfTerms(n2_new,dim=3)
            B_new = numpy.zeros((2,nn2_new,ntilt))

            for itilt in range(ntilt):
                proj_map = util.PolynomialTransformation3D(self.resid.n2,B[:,:,itilt])
                proj_map = proj_map.compose(flat_to_warp_tf)
                B_new[:,:,itilt] = proj_map.coeffs()[:2,:]

            self.resid.extendProjectionMapOrder( n2_new, full_reset=False )

            self.resid.n2 = n2_new
            self.resid.nn2 = nn2_new

            self.b = B_new.ravel().tolist()

        else:   # Relax the projection map

            self.resid.a = self.a

            n2_new = self.flatten_order_eff*self.resid.n2

            if n2_new>self.n2:
                self.resid.extendProjectionMapOrder( n2_new, full_reset=False )

            self.resid.b = self.resid.relaxProjection()

            self.b = self.resid.b

        
       # if self.doPlot:
        if self.showPlot3D:
            from enthought.mayavi import mlab
            mlab.show()



    def reFrameSeries(self):
        '''Re-evaluate the natural reconstruction frame (XY grid boundaries) for each series.
        Calulates by back-projecting the four corners of the reference tilt at Z=0
        '''
        
        log.info("Evaluate XY output frame:")
        
        minimums = []
        maximums = []
        
        ntilt = self.numberOfTilts
        nn2 = self.resid.nn2
        
        b = numpy.resize(self.b,(2,nn2,ntilt))
        
        for index,series_name in enumerate(self.seriesBasenames):
            
            series = self.project.getSeries(series_name)
            
            (nx_src,ny_src) = series.dimensions()
            
            indexOfReferenceExposure = self.referenceExposures[index]
            
            b_ = b[:,:3,indexOfReferenceExposure]    # We disregard the Z coordinate (Z=0)
            
            tf = util.AffineTransformation2D(b_).inv()
            
            if series.input=='st':
                T = series.prealignTranslations[series.indexOfReferenceExposure(),:]
                translation = util.AffineTransformation2D( numpy.column_stack((T,numpy.eye(2))))
                tf = tf.compose(translation) 
            
            corners_src = numpy.array([[0,0],[0,ny_src-1],[nx_src-1,0],[nx_src-1,ny_src-1]])
            corners_img = tf.forward(corners_src)
                                    
            series_minimums = numpy.min(corners_img,axis=0).astype('int')
            series_maximums = numpy.max(corners_img,axis=0).astype('int')
            
            log.info( "Series #%i: (xmin,xmax)=(%i,%i)  (ymin,ymax)=(%i,%i)" \
                      %(index,series_minimums[0],series_maximums[0],series_minimums[1],series_maximums[1]) )
            
            series.xmin = series_minimums[0]
            series.xmax = series_maximums[0]
            series.ymin = series_minimums[1]
            series.ymax = series_maximums[1]
            series.zmin = self.project.reconstruction.origin.z
            series.zmax = self.project.reconstruction.end.z
            
            minimums.append(series_minimums)
            maximums.append(series_maximums)
            
            for icorner in range(4):
                log.debug("(%.2f,%.2f) -> (%.2f,%.2f)" \
                          %(corners_src[icorner,0],corners_src[icorner,1],corners_img[icorner,0],corners_img[icorner,1]))

        minimums = numpy.min(numpy.row_stack(minimums),axis=0)
        maximums = numpy.max(numpy.row_stack(maximums),axis=0)

        self.project.nx = maximums[0] - minimums[0] + 1
        self.project.ny = maximums[1] - minimums[1] + 1
        
        self.project.reconstruction.setOrigin(minimums[0], minimums[1], None)
        self.project.reconstruction.setEnd(maximums[0], maximums[1], None)



    def evaluateRemapTransform( self ):
        '''Evaluate the 2D remap transformation to be used during the filtering.'''
        
        # This routine will only use the point-like objects to evaluate 
        # the remap
        
        log.info("Evaluate remap function used during the filtering!")
        
        log.warning("REMAP ST OR PREALI")
        
        # Gather the fiducial informations in 2 and 3D
                
        #XYZ,obj_mask = self.resid.extract3DPointMarkers()
        XYZ,obj_mask = self.extract3DPointMarkers()

#        pointTracks = numpy.squeeze(self.tracks[:,self.points,:,:1])   # tracks[itilt,ipoint,axis,order]
#        pointMask = self.mask[:,self.points]    # mask[itilt,ipoint]
        pointTracks = numpy.squeeze(self.tracks[:,obj_mask,:,:1])   # tracks[itilt,ipoint,axis,order]
        pointMask = self.mask[:,obj_mask]    # mask[itilt,ipoint]
                

        # Do the remap for each tilt series that is enabled
        
        for index,series_name in enumerate(self.seriesBasenames):
            
            series = self.project.getSeries(series_name)
                        
            series_rg = self.seriesRange[index]
            
            tracks = pointTracks[series_rg,:,:]
            mask = pointMask[series_rg]                
            
            # Take into account for the relative rotation of the sample between each tilt series
            
            tf = series.getRelativeReferenceTransform()
            XYZ_tf = tf.forward(XYZ)
                        
            # Gather the geometrical tilt series information
            
            M0 = numpy.asarray(self.M0[index])
            u_t = numpy.asarray(self.u_t[index])
            angles = self.angles[1,series_rg] # Replace b the truly calculated angles
            
            log.debug("\nSeries %s:" %series_name)
            log.debug("   %.30s ..." %str(numpy.asarray(series_rg)))
            log.debug("   %.30s ..." %numpy.asarray(angles))
            log.debug("   %s" %tf)
        
            # Call the remap evaluation procedure - Results will be stored in
            # a file created from the basename and model in the filter_directory
            
            nx,ny = series.dimensions()
            center = numpy.array([nx/2.0,ny/2.0,0.0])
            
            txbr.filter.remap( tracks, mask, XYZ_tf, M0, u_t, angles, \
                               self.remap_order, self.remap_magnification, \
                               doPlot=self.extraPlot or True, directory=self.filter_directory,
                               basename=series_name, model=self.model(), center=center )


    def evaluateAlignmentTransform( self  ):
        '''In this routine, the 2D alignment transform is calculated. This 
        transformation allows to generate images similar to the .ali files
        of etomo. The source file can be either the st or preali files
        depending on the input used.'''
        
        log.info('Calculate alignment transform...')
                
        nn2 = self.resid.nn2
        ntilt = self.numberOfTilts

        # We use the ideal projection map as a reference to align the micrographs between
        # themselves

        if self.mosaic:
            set = 'global'  # single global rotation axis (image shift mosaic case)
        else:
            set = 'local'  # ex-centered rotation axis (mostly stage shift mosaic case)

     #   print set

        angles_ideal, P_ideal = self.initializeProjectionMap( mode='ideal', set=set )

        b_ortho = numpy.zeros((2,nn2,ntilt))
        b_ortho[0,:4,:] = P_ideal[::2,:4].T
        b_ortho[1,:4,:] = P_ideal[1::2,:4].T
        
        # Relate to the optimal projection map
        
        b = numpy.resize(self.b,(2,nn2,ntilt))
        
        R = numpy.zeros((2,3,ntilt))
        
        for itilt in range(ntilt):
            
            M = numpy.zeros((4,3))
            
            M[0,0] = 1.0
            M[0,1:3] = b[:,0,itilt]
            M[1,1:3] = b[:,1,itilt]
            M[2,1:3] = b[:,2,itilt]
            M[3,1:3] = b[:,3,itilt]
            
            Minv = numpy.linalg.pinv(M)    
            
            b_ortho_x =  b_ortho[0,:4,itilt]
            b_ortho_y =  b_ortho[1,:4,itilt]
            
            R[0,:,itilt] = numpy.dot(Minv,b_ortho_x)
            R[1,:,itilt] = numpy.dot(Minv,b_ortho_y)
            
        log.info('Save the alignment tranforms to create an align stack...')

        nseries = self.project.numberOfEnabledSeries()

        for iseries in range(nseries):

            rg = self.seriesRange[iseries]
            
            filename = os.path.join(self.align_directory,'%s.tf2ali' %(self.seriesBasenames[iseries]))
             
            f = open(filename,'w')
            
            for itilt in rg:
                f.write('%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n' %(R[0,1,itilt],R[0,2,itilt],R[1,1,itilt],R[1,2,itilt],R[0,0,itilt],R[1,0,itilt]))
            
            f.close()


    def initializeStructureSet( self, objects=None ):
        '''Initialize the Structure Set object'''

        log.debug("Initialize the Structure Set object.")
        
        npatch = self.numberOfPatches
        ntilt = self.numberOfTilts

        nn1 = self.resid.nn1
        nn4 = self.resid.nn4

        B = (self.resid.n2,numpy.array(self.b))

        xmax = self.project.nx
        ymax = self.project.ny

        self.structureSet = modl.StructureSet( B, xmax, ymax, directory=self.align_directory, \
                                               basename=self.project.basename )

        for ipatch in range(self.numberOfPatches):

            if objects!=None and objects.count(ipatch)==0:
                continue

            A = (self.resid.n1,numpy.array(self.a[ipatch:len(self.a):npatch]))
            C = (self.resid.n3,numpy.array(self.c[ipatch:len(self.c):npatch]))
            D = (self.resid.n4,numpy.array(self.d[ipatch:len(self.d):npatch]))

            pt = ipatch in self.points
            ln = ipatch in self.lines

            structure = modl.Structure( A, C, D, isAPoint=pt, isALine=ln )
            self.structureSet.addNewStructure(structure)

            if not ipatch in self.points:    # Use projective duality to initialize the structure
                log.info('Initialize structure %i/%i' %(ipatch+1,self.numberOfPatches))
                parameters = self.glomap[ipatch]
                keywords = { 'directory': self.contour_directory, 'basename': self.project.basename }
                structure.initializeSurfaceMapping( *parameters, **keywords )
                self.a[ipatch:len(self.a):npatch] = structure.a[:]

        if self.showPlot3D and not self.pure_bead_markers:
            self.structureSet.plot_structures(direct=True) # Plot Estimates for the 3D patches
            self.structureSet.plot_structures()

        file = os.path.join( self.align_directory, self.project.basename + '.init.mod')
        
        nx = self.project.nx
        ny = self.project.ny
        nz = int(min(nx,ny)/2)
        
        sx = self.project.reconstruction.scale.x
        sy = self.project.reconstruction.scale.y
        sz = self.project.reconstruction.scale.z        
        
        self.structureSet.saveModel( self.project.nx, self.project.ny, \
                                     int(min(self.project.nx,self.project.ny)/2), \
                                     scalex=sx, scaley=sy, scalez=sz, 
                                     show=False, filename=file, direct=True)


    def saveStructureSet(self):

        nn1 = self.resid.nn1
        nn4 = self.resid.nn4

        B = (self.resid.n2,numpy.array(self.b))

        self.structureSet.nb,self.structureSet.b = B

        npatch = self.numberOfPatches
        ntilt = self.numberOfTilts

        for ipatch in range(self.numberOfPatches):

            A = (self.resid.n1,numpy.array(self.a[ipatch:len(self.a):npatch]))
            C = (self.resid.n3,numpy.array(self.c[ipatch:len(self.c):npatch]))
            D = (self.resid.n4,numpy.array(self.d[ipatch:len(self.d):npatch]))

            structure = self.structureSet.structures[ipatch]

            structure.update(A,C,D)

#        nx = self.project.nx
#        ny = self.project.ny
#        nz = int(min(nx,ny)/2)
        
        nx = self.project.reconstruction.end.x - self.project.reconstruction.origin.x
        ny = self.project.reconstruction.end.y - self.project.reconstruction.origin.y
        nz = self.project.reconstruction.end.z - self.project.reconstruction.origin.z
        
        log.warning("Check the Z-offset for the model")
        
        offset = self.project.reconstruction.getOffset()
        
        xoffset = offset.x
        yoffset = offset.y
        zoffset = offset.z

        xoffset = 0.0
        yoffset = 0.0
        zoffset = 0.0
        
        sx = self.project.reconstruction.scale.x
        sy = self.project.reconstruction.scale.y
        sz = self.project.reconstruction.scale.z        

        self.structureSet.saveModel( nx, ny, nz, xoffset, yoffset, zoffset, \
                                     scalex=sx, scaley=sy, scalez=sz, show=self.showPlot3D )

        log.info(self.structureSet.info())


    def saveProjectionMapToProject( self, n2=None, nterm=None, b2proj=None, extension=None ):
        """
        Copy the projection map within the project, and generate the configuration files.
        Projection maps stored under their general expanded expressions.

        :param n2: The order of the projection maps.
        :param nterm: The number of terms in a projection polynomial function.
        :param b2proj: The projection map coefficients.
        :param extension: An extension that will be added to the configuration file names.
        """

        log.info("Save the new projection maps!")

        ntilt = self.numberOfTilts

        if n2==None or nterm==None or b2proj==None:
            n2 = self.resid.n2
            nterm = self.resid.nn2
            b2proj = self.b

        for index,basename in enumerate(self.seriesBasenames):

            series = self.project.getSeries(basename)
            series.projection.approximationOrder(approximationOrder=n2)

            rg = self.seriesRange[index]

            for iproj in range(nterm):
                for itilt in rg:
                    iexp = (itilt-rg[0])*self.decimate
                    series.projection.x_coefficients[iexp][iproj] = b2proj[iproj*ntilt+itilt]
                    series.projection.y_coefficients[iexp][iproj] = b2proj[nterm*ntilt+iproj*ntilt+itilt]
                series.projection.x_coefficients[:,iproj] = __re_sample_coefficients__(series.projection.x_coefficients[:,iproj],self.decimate)
                series.projection.y_coefficients[:,iproj] = __re_sample_coefficients__(series.projection.y_coefficients[:,iproj],self.decimate)

            for itilt in rg:
                iexp = (itilt-rg[0])*self.decimate
                series.projection.scaling_coefficients[iexp][1:] = 0.0
                for i in range(4):
                    series.projection.scaling_coefficients[iexp][i] = self.b_scaling[i*ntilt+itilt]

            for iproj in range(4):
                series.projection.scaling_coefficients[:,iproj] = __re_sample_coefficients__(series.projection.scaling_coefficients[:,iproj],self.decimate)

        for index,basename in enumerate(self.seriesBasenames):
            if series.input=='st':
                log.info('Input from st file: take a prealignment into account!')
                series.projection.x_coefficients[:,0] -= series.prealignTranslations[:,0]
                series.projection.y_coefficients[:,0] -= series.prealignTranslations[:,1]
                    
        self.project.reconstruction.alignment_model = self.__model_representation__( n2=n2 )
        
        self.project.saveAll( extension=extension )


    def loadParameters(self,file):
        """
        Load the a,b,c,d parameters from a file.
        """

        try:
            f = open(file)
            try:
                a = f.readline()
                b = f.readline()
                b_scaling = f.readline()
                c = f.readline()
                d = f.readline()
            finally:
                f.close()
        except IOError:
            pass

        a = [float(item) for item in a.split('\t')]
        b = [float(item) for item in b.split('\t')]
        b_scaling = [float(item) for item in b_scaling.split('\t')]
        c = [float(item) for item in c.split('\t')]
        d = [float(item) for item in d.split('\t')]

        return [a,b,b_scaling,c,d]


    def saveParameters( self, file ):
        """
        Write the a,b,c,d coefficient parameters into a file.

        :param file: The path of the file where to store the coefficients.
        """

        f = open(file,'w')
        try:
            f.write('\t'.join(str(i) for i in self.resid.a) + '\n')
            f.write('\t'.join(str(i) for i in self.resid.b) + '\n')
            f.write('\t'.join(str(i) for i in self.resid.b_scaling) + '\n')
            f.write('\t'.join(str(i) for i in self.resid.c) + '\n')
            f.write('\t'.join(str(i) for i in self.resid.d) + '\n')
        finally:
            f.close()


    def info(self):
        """This function calculates and displays the different error terms (reprojection
        error, tangency error) in the cost function. It can be easily evaluated with respect
        to the tilt, patch and so on...
        """
        
        log.info('Residual Information')

        t_start = time.time()

        ntilt = self.numberOfTilts
        npatch = self.numberOfPatches

        e = self.resid.value()

        if e==None:  # Cannot do the adequate plots
            log.info('Skip info routine!')
            return

        e1,e2 = e
        res1,res2 = numpy.sqrt(e1),numpy.sqrt(e2)
        res = res1 + res2
        
        E1 = e1.sum()
        E2 = e2.sum()
        
        d = res1.sum()/self.mask.sum()
        #std = numpy.std(e1)
        std = numpy.sqrt(e1.sum()/self.mask.sum())

        for itilt in range(ntilt):
            for ipatch in range(npatch):
                e1_ = e1[itilt,ipatch]
                e2_ = e2[itilt,ipatch]
                log.debug('itilt=%-3i ipatch=%-3i  E=%6.2e  E1=%6.2e  E2=%6.2e' %(itilt,ipatch,e1_+e2_,e1_,e2_))
                #print 'itilt=%-3i ipatch=%-3i  E=%6.2e  E1=%6.2e  E2=%6.2e' %(itilt,ipatch,e1_+e2_,e1_,e2_)

        log.info('                     Total Error:  %6.2e = %6.2e + %6.2e' %(E1+E2,E1,E2))
        log.info('           Error by tilt & patch:  %6.2e = %6.2e + %6.2e' %((E1+E2)/ntilt/npatch,E1/ntilt/npatch,E2/ntilt/npatch))
        log.info('    Error by tilt & patch & axis:  %6.2e = %6.2e + %6.2e' %((E1+E2)/ntilt/npatch/2.0,E1/ntilt/npatch/2.0,E2/ntilt/npatch/2.0))
        log.info('Mean Reprojection Distance Error:  %6.2e' %d)
        log.info('        Standard Deviation Error:  %6.2e\n' %std)

        for itilt in range(ntilt):
            E1byTilt = e1.sum(axis=1)
            E2byTilt = e2.sum(axis=1)
            EbyTilt = E1byTilt + E2byTilt
            log.info('  Tilt: %3i/%-3i    E1=%-10.2f  E2=%-10.2f  All Patches: %s' \
                     %(itilt,ntilt,E1byTilt[itilt],E2byTilt[itilt],numpy.all(self.mask[itilt,:])))

        for ipatch in range(npatch):
            E1byPatch = e1.sum(axis=0)
            E2byPatch = e2.sum(axis=0)
            EbyPatch = E1byPatch + E2byPatch
            log.info(' Track: %3i/%-3i    E1=%-10.2f  E2=%-10.2f  All Tilts: %s' \
                     %(ipatch,npatch,E1byPatch[ipatch],E2byPatch[ipatch],numpy.all(self.mask[:,ipatch])))

        t_end = time.time()

        log.debug('Time to calculate direct Error: %f' %(t_end-t_start))
        
        # Do some plottings

        errorByTilt = True and self.doPlot
        errorByPatch = True and self.doPlot
        errorByTiltAndPatch = True and self.doPlot
                   
        if errorByTilt: 
            txbr.utilities.plotErrorByTilt( E1byTilt, E2byTilt, EbyTilt, \
                                            directory=self.align_directory, basename=self.project.basename, \
                                            model = self.model() )
        
        if errorByPatch: 
            txbr.utilities.plotErrorByPatch( E1byPatch, E2byPatch, EbyPatch, \
                                             directory=self.align_directory, basename=self.project.basename, \
                                             model = self.model() )
 
        if errorByTiltAndPatch: 
            txbr.utilities.plotErrorByTiltAndPatch( res1, res2, res, \
                                             directory=self.align_directory, basename=self.project.basename, \
                                             model = self.model() )
            
        # dump info in a file
        
        dump_info_in_file = True
        
        if dump_info_in_file:

            current_residual_filename = os.path.join( self.align_directory, "residual_%s.txt" %self.model())
            filename = os.path.join( self.align_directory, "residual.txt")
            
            if os.path.lexists(filename): os.unlink(filename)
            
            os.symlink(os.path.basename(current_residual_filename),filename)

            info_file = open(current_residual_filename,"w")
            
            info_file.write('Alignment: %s\n' %self.__model_representation__())
            info_file.write('\n')
            
            info_file.write('Total Error:  %6.2e = %6.2e + %6.2e\n' %(E1+E2,E1,E2))
            info_file.write('Error by tilt & patch:  %6.2e = %6.2e + %6.2e\n' %((E1+E2)/ntilt/npatch,E1/ntilt/npatch,E2/ntilt/npatch))
            info_file.write('Error by tilt & patch & axis:  %6.2e = %6.2e + %6.2e\n' %((E1+E2)/ntilt/npatch/2.0,E1/ntilt/npatch/2.0,E2/ntilt/npatch/2.0))
            info_file.write('Mean Reprojection Distance Error:  %6.2e\n' %d)
            info_file.write('Standard Deviation Error:  %6.2e\n' %std)
                    
            info_file.write('\n')
            info_file.write('Error by track:\n')
            
            for ipatch in range(npatch):
                E1byPatch = e1.sum(axis=0)
                E2byPatch = e2.sum(axis=0)
                EbyPatch = E1byPatch + E2byPatch
                info_file.write('Track: %3i/%-3i    E1=%-10.2f  E2=%-10.2f  All Tilts: %s\n' \
                         %(ipatch,npatch,E1byPatch[ipatch],E2byPatch[ipatch],numpy.all(self.mask[:,ipatch])))
            
            info_file.write('\n')
            info_file.write('Error by tilt:\n')

            for itilt in range(ntilt):
                E1byTilt = e1.sum(axis=1)
                E2byTilt = e2.sum(axis=1)
                EbyTilt = E1byTilt + E2byTilt
                info_file.write('Tilt: %3i/%-3i    E1=%-10.2f  E2=%-10.2f  All Patches: %s\n' \
                         %(itilt,ntilt,E1byTilt[itilt],E2byTilt[itilt],numpy.all(self.mask[itilt,:])))

            info_file.write('\n')
            info_file.write('Error by tilt and track:\n')
            
            for itilt in range(ntilt):
                for ipatch in range(npatch):
                    e1_ = e1[itilt,ipatch]
                    e2_ = e2[itilt,ipatch]
                    info_file.write('Tilt: %3i/%-3i Track: %3i/%-3i  E=%6.2e  E1=%6.2e  E2=%6.2e\n' %(itilt,ntilt,ipatch,npatch,e1_+e2_,e1_,e2_))

            info_file.close()
                        
            # Sorted Error
            
            indices = list(numpy.ndenumerate(res))

            current_residual_filename = os.path.join( self.align_directory, "residual_sorted_%s.txt" %self.model())
            filename = os.path.join( self.align_directory, "residual_sorted.txt")
            
            if os.path.lexists(filename): os.unlink(filename)
            
            os.symlink(os.path.basename(current_residual_filename),filename)

            sorted_info_file = open(current_residual_filename,"w")
            
            arg = numpy.argsort(res,axis=None)
            
            sorted_info_file.write("Tilts and patches are indexed from 0...\n")
            for index in arg[::-1]:
                (itilt,ipatch), err = tuple(indices[index])
                sorted_info_file.write("Tilt #%-5i    Patch #%-5i    ->    %.2f\n" %(itilt, ipatch, err))
            
            sorted_info_file.close()
            
        # Log info into a file
            
        log_info_in_file = True
        
        if log_info_in_file:

            filename = os.path.join( self.align_directory, "txbr.align.bundle.log" )
            
            log_file = open(filename,"a")
            
            log_file.write(time.strftime('-------------------- %a, %d %b %Y %H:%M:%S --------------------\n',time.gmtime()))
            log_file.write('Alignment: %s\n' %self.__model_representation__())
            log_file.write('Series: %s\n' %(", ".join(self.seriesBasenames)))
            log_file.write('Order: (n1,n2,n3,n4)=(%i,%i,%i,%i)\n' %(self.n1,self.n2,self.n3,self.n4))
            log_file.write('Total Error:  %6.2e = %6.2e + %6.2e\n' %(E1+E2,E1,E2))
            log_file.write('Error by tilt & patch:  %6.2e = %6.2e + %6.2e\n' %((E1+E2)/ntilt/npatch,E1/ntilt/npatch,E2/ntilt/npatch))
            log_file.write('Error by tilt & patch & axis:  %6.2e = %6.2e + %6.2e\n' %((E1+E2)/ntilt/npatch/2.0,E1/ntilt/npatch/2.0,E2/ntilt/npatch/2.0))
            log_file.write('Mean Reprojection Distance Error:  %6.2e\n' %d)
            log_file.write('Standard Deviation Error:  %6.2e\n' %std)

            log_file.close()
       

    def generateAlignLogFile(self,threshold_error=3.0):
        '''Geneate an align log file for the point-like markers of each series'''
        
        if len(self.lines)!=0 or len(self.surfaces)!=0: return
        
        log.info('Generate the alignment log file!')

        # Build the data array
        # Array containing: itilt   ipatch  x   y   xp  yp  (xp-x)  (yp-y)  std

        data = numpy.zeros((self.numberOfTilts,len(self.points),7))

        tilts = numpy.repeat(numpy.arange(self.numberOfTilts),len(self.points))
        tilts = numpy.resize(tilts,(self.numberOfTilts,len(self.points)))

        patches = numpy.arange(len(self.points))    # would be broadcasted

        tracks_ = numpy.rollaxis(self.tracks[:,:,:2,0],2,0)*self.mask
        proj_ = self.resid.directProject()*self.mask

        data[:,:,0] = tilts
        data[:,:,1] = patches
        data[:,:,2] = tracks_[0,:,:]
        data[:,:,3] = tracks_[1,:,:]
        data[:,:,4] = proj_[0,:,:] - tracks_[0,:,:]
        data[:,:,5] = proj_[1,:,:] - tracks_[1,:,:]
        data[:,:,6] = numpy.sqrt(data[:,:,4]**2 + data[:,:,5]**2)

        # Generate the log file

        for iseries,basename in enumerate(self.seriesBasenames):

            error_file = "align." + basename + ".log"
            error_file = os.path.join(self.align_directory,error_file)

            log.info('... Generate %s' %error_file)
            
            f = open(error_file,'w')
            
            f.write('         Projection points with large residuals\n')
            f.write(' obj  cont  view   index coordinates    residuals    # of\n')
            f.write('   #     #     #      X         Y        X      Y    S.D.\n')
            
            tilt_range = self.seriesRange[iseries]
            ntilt = self.tlts[iseries]
            
            data_series = numpy.resize(data[tilt_range],(ntilt*len(self.points),7))
            
            indices = numpy.argsort(data_series[:,6])
            data_series = data_series[indices[::-1]]
            
            std_dev = numpy.std(data_series[:,6])
            
            for s in data_series:
                if s[6]<threshold_error*std_dev: 
                    continue
                itilt,ipatch,X,Y,resX,resY = s[:6]
                values = {'object':1}
                values['contour'] = ipatch + 1
                values['view'] = itilt + 1 - tilt_range[0]
                values['X'] = X
                values['Y'] = Y
                values['resX'] = resX
                values['resY'] = resY
                values['SD'] = s[6]/std_dev
                f.write('%(object)4.0f %(contour)5.0f %(view)5.0f %(X)9.2f %(Y)9.2f %(resX)6.2f %(resY)6.2f %(SD)6.2f\n' %values)
            
            f.close()
 

    def plot_individual_contour_mappings(self):

        log.info('Plot individual contour mapping')

        # First build the array containing the experimental points

        xy_exp = []

        for mod in self.models:
            for object in mod.objects:
                for contour in object.contours:
                    xy_exp.append(contour.points[:,:2])

        xy_exp = numpy.row_stack(xy_exp)
        (x_exp,y_exp) =  numpy.array_split(xy_exp,2,axis=1)

        x_exp = x_exp.ravel()
        y_exp = y_exp.ravel()

        # Second build the array containing the experimental points

        xy_ct_pts,xy_ct = [],[]

        n3 = self.resid.n3
        nn3 = self.resid.nn3
        orders = [n3]

        t_pt = numpy.zeros((1,1))

        t_surf = numpy.linspace(0.0,1.0,num=20)
        t_surf.resize((t_surf.size,1))

        for mod in self.models:
            index0 = mod.tiltOffset
            for object in mod.objects:
                iobject = object.indexOfObject
                if object.isAPoint():
                    t = t_pt
                else:
                    t = t_surf
                for contour in object.contours:
                    #itilt = index0 + contour.indexOfContour
                    z = int(round(contour.zExtensionSet().pop(),0))
                    itilt = index0 + z
                    coeff_x = self.tracks[itilt,iobject,0,:nn3]
                    coeff_y = self.tracks[itilt,iobject,1,:nn3]
                    pc_x = util.Poly(coeff_x,orders)
                    pc_y = util.Poly(coeff_y,orders)
                    if object.isAPoint():
                        xy_ct_pts.append(pc_x.eval(t))
                        xy_ct_pts.append(pc_y.eval(t))
                    else:
                        xy_ct.append(pc_x.eval(t))
                        xy_ct.append(pc_y.eval(t))

        if len(xy_ct)!=0:
            xc = numpy.column_stack(xy_ct[::2])
            yc = numpy.column_stack(xy_ct[1::2])
            xy_ct = numpy.row_stack((xc,yc))
            (x_ct,y_ct) = numpy.array_split(xy_ct,2,axis=0)
        else:
            x_ct,y_ct = numpy.zeros((0)),numpy.zeros((0))

        if len(xy_ct_pts)!=0:
            xc_pts = numpy.column_stack(xy_ct_pts[::2])
            yc_pts = numpy.column_stack(xy_ct_pts[1::2])
            xy_ct_pts = numpy.row_stack((xc_pts,yc_pts))
            (x_ct_pts,y_ct_pts) = numpy.array_split(xy_ct_pts,2,axis=0)
        else:
            (x_ct_pts,y_ct_pts) = numpy.zeros((0)),numpy.zeros((0))

        # do the plotting

        nx,ny = self.project.nx,self.project.ny

        f = pylab.figure()

        pylab.scatter(x_exp,y_exp,c='green')
        pylab.plot(x_ct_pts,y_ct_pts,'g:')    # contours
        pylab.plot(x_ct,y_ct,'g:')    # Point-like contours

        pylab.xlim(0,nx)
        pylab.ylim(0,ny)

        #pylab.title(r'Individual Contour mapping')
        pylab.xlabel(r'$x$')
        pylab.ylabel(r'$y$')
        
        pylab.title(r'Individual Contours')


    def plot_projections(self):
        '''This routine plots the projection coefficients as a function of the tilt
        index'''

        ntilt = self.numberOfTilts
        nterm = self.resid.nn2
  
        P_ = numpy.resize(self.b,(2,nterm,ntilt))
        
        P_id = numpy.resize(self.P,(ntilt,2,4))
        P_id = numpy.rollaxis(P_id,0,3) 
    
        txbr.utilities.plot_projections( P_, P_id, directory=self.align_directory, \
                                         basename=self.project.basename, model = self.model() )
        
        
    def plot_angles(self, angles=None):
        '''If the residual object makes an orthogonal asumption, then the angles intervening
        in the linear part of the projection map are plotted.'''
        
        try:
        
            ntilt = self.numberOfTilts
            nterm = self.resid.nn2
        
            angles_id = __evaluate_orthogonal_from_projection__(self.P)
        
            if angles==None:
                angles = self.resid.getOrthogonalParameters()
            else:
                angles = numpy.asarray(angles)
                angles = numpy.resize(angles,angles_id.shape)
                        
            txbr.utilities.plot_angles(angles, angles_id)
            
        except NotImplementedError:
            
            pass
        
        
    def plot_magnifications( self ):
        '''If the residual object makes an orthogonal asumption, then the angles intervening
        in the linear part of the projection map are plotted.'''
        
        try:
        
            txbr.utilities.plot_magnifications(self.resid.getMagnificationParameters())
            
        except NotImplementedError:
            
            pass


    def plot_contour_reprojection( self, itilt=None, ipatch=None, t0=0.0, t1=1.0, limits=None ):
        """This function displays the reprojection of a patch for a given tilt. 
        If no tilt is specified, the reprojection for the available tilts will be displayed,
        but plots will be separated in as many enabled series. 
        If no patch is specified, the reprojection for the patches is displayed.
        Variable <itilt>: specifies the tilt index to be displayed
        Variable <ipatch>: specifies the patch be displayed
        """

        log.info('Plot Contour Reprojections.')

        try:
            self.resid.Pvars
            self.resid.P1
            self.resid.P2
        except AttributeError:
            log.debug("Problem with the resid object: no Pvars, P1 or P2")
            return

        t_start = time.time()
        
        # Determine what are the series that should be involved for this plot
        
        if itilt==None:
            series_ = [ self.project.getSeries(basename) for basename in self.seriesBasenames]
            models_ = self.models
            seriesRange_ = self.seriesRange
        else:
            series_ = []
            models_ = []
            seriesRange_ = []
            for index,basename in enumerate(self.seriesBasenames):
                if itilt in self.seriesRange(index):
                    series_.append(self.project.getSeries(basename))
                    models_.append(self.models[index])
                    seriesRange_.append([itilt])
        
        # Determine what are the object (patches) that should be involved for this plot
                
        if ipatch==None:
            patches = range(self.numberOfPatches)
        else:
            patches = [ ipatch ]
            
        nobj = len(patches)
                    
        # Get the marker traces . 
        # Stored in a list of list of list: nseries->nobj->npy_arr(n,2) with n less than ntilt*ndata
        traces = [[ model.objects[iobj_].points() for iobj_ in patches ] for model in models_ ]
        
        # Evaluate the re-projection
        # Stored in a list of list of list: nseries->nobj->ntilt->npy_arr(ndata,2)
        # ndata may depend on the object type
        if len(self.surfaces)==0 and len(self.lines)==0:
            reproj = self.resid.directProject()
            reproj = numpy.swapaxes(reproj,0,2)
            reproj = numpy.resize(reproj,(self.numberOfPatches,1,self.numberOfTilts,2))
            reproj = numpy.array([ reproj[:,:,rg,:] for rg in seriesRange_ ])
            reproj = reproj[:,patches,:,:]
        else:
            reproj = [[[ numpy.column_stack(self.resid.project(itilt_,iobj_)) for itilt_ in rg ] for iobj_ in patches ] for s,rg in zip(series_,seriesRange_) ]


        # Do the figure
        
        for index,s in enumerate(series_):
            
            if limits==None: limits = s.getSourceFrame(full=False).ravel()
            
            marks = numpy.row_stack(traces[index])

            pylab.figure()

            pylab.scatter(marks[:,0],marks[:,1],c = 'green',label='Tracks')
            for eval_reproj in reproj[index]:
                rep = numpy.asarray(eval_reproj)
                pylab.plot(rep[:,:,0], rep[:,:,1], 'r:1', label='Reprojected Data')
        
            title = s.basename.replace("_","\_")
                
            pylab.axis(limits)
            pylab.xlabel(r'x')
            pylab.ylabel(r'y')
            pylab.title(title)
            pylab.legend(('Reprojected Data',))
            
            if itilt!=None and ipatch!=None:
                extension = "_tilt-%i_patch-%i" %(itilt+1,ipatch+1)
                pylab.figtext(0.15,0.85,'Patch \#%i/%i  Tilt \#%i/%i' %(ipatch+1,npatch,itilt+1,ntilt),fontsize=14)
            elif itilt!=None:
                extension = "_tilt-%i" %(itilt+1)
                pylab.figtext(0.15,0.85,'Tilt \#%i/%i' %(itilt+1,self.numberOfTilts-1),fontsize=14)
            elif ipatch!=None:
                extension = "_patch-%i" %(ipatch+1)
                pylab.figtext(0.15,0.85,'Patch \#%i/%i' %(ipatch+1,self.numberOfPatches),fontsize=14)
            else:
                extension = ""
    
            current_file = "%s%s_%s_reproj.png" %( s.basename, extension, self.model() )
            current_file = os.path.join( self.project.align_directory, current_file )
    
            pylab.savefig( current_file, format='png' )
            
            file = "%s%s_reproj.png" %( s.basename, extension )
            file = os.path.join( self.project.align_directory, file )
            
            if os.path.lexists(file): os.unlink(file)
            
            os.symlink(os.path.basename(current_file),file)
            
        t_end = time.time()
    
        log.debug('Time to do the reprojection plot: %f' %(t_end-t_start))
            
        t_end = time.time()
    
        log.debug('Time to do the reprojection plot: %f' %(t_end-t_start))


    def plot_Initial_Estimates_Reprojection(self):
        """Reproject the tangent points of the surface with the electron beam
        obtained with the Dual Space Approach.
        """

        xy = []

        for mod in self.models:
            index0 = mod.tiltOffset
            for object in mod.objects:
                iobject = object.indexOfObject
                for contour in object.contours:
                    for point in contour.points:
                        indexOfTilt = index0 + int(round(point[2],0))
                        xy.append(point[:2])

        (xdata,ydata) = numpy.array_split(numpy.array(xy),2,axis=1)

        xdata = numpy.resize(xdata,(xdata.size))
        ydata = numpy.resize(ydata,(ydata.size))

        x,y = [],[]

        for structure in self.structureSet.structures:
            try:
                indicesOfTilts = structure.indicesOfTilts
                XYZ = structure.points_OC
            except AttributeError:
                continue
            (ntilt,ndata,dim) = XYZ.shape
            for itilt in range(ntilt):
                indexOfTilt = indicesOfTilts[itilt]
                x.append(self.P[2*indexOfTilt,0] + numpy.tensordot(XYZ[itilt,:,:],self.P[2*indexOfTilt,1:4],axes=([1],[0])))
                y.append(self.P[2*indexOfTilt+1,0] + numpy.tensordot(XYZ[itilt,:,:],self.P[2*indexOfTilt+1,1:4],axes=([1],[0])))

        # Do the plot

        pylab.figure()

        pylab.scatter(xdata, ydata)
        for itilt in range(len(x)):
            pylab.plot(x[itilt],y[itilt],'r:')

        pylab.title(r'Reprojection of Tangent Points obtained from Dual Space Approach',fontsize=16,color='r')
        pylab.xlabel(r'x',fontsize=18)
        pylab.ylabel(r'y',fontsize=18)

        nx,ny = self.project.nx,self.project.ny

        pylab.xlim(0,nx)
        pylab.ylim(0,ny)

        file = self.project.basename + '_reproject_DSA.png'
        file = os.path.join(self.project.work_directory,align_dir,file)
        pylab.savefig(file,format='png')