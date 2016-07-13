import sys
import os
import numpy
import scipy.signal
import mrc
import txbr.utilities
import txbr.prj

from txbr.onthefly import QuickReconstruction

class CrossReconstruction(QuickReconstruction):

    def bckprj( self, limits=None, crossvol=None, weight=None, output=None, showVol=False ):
        '''Back-project the slices into the final volume'''

        if output==None: output = "%s_vol.mrc" %(self.stack.basename)

        if len(self.angles2process)==0: return

        if limits!=None:
            self.Zmin = int(limits[0])
            self.nz = int(limits[1]-limits[0]+1)

        mrcVolume = self.createVolume( output=output, from_scratch=len(self.processed)==0 )

        volumes = []

        for iz in range(self.nz):
            slice = mrcVolume.getZSliceAt(iz)
            volumes.append(slice)

        src_file = "%s.SL" %( self.stack.basename )
        fsrc = mrc.MRCFile(src_file)

        if weight!=None: fweight = mrc.MRCFile(weight)
        if crossvol!=None:
            fcross = mrc.MRCFile(crossvol)

        for index,angle in enumerate(self.angles2process):

            sys.stdout.write( "Angle: %i: Z-> " %int(angle) )

            u = fsrc.getZSliceAt(index)
            u = u - numpy.mean(u)
            print "u: mean: %.2f std: %.2f" %(numpy.mean(u),numpy.std(u))
            if weight!=None:
                weight = fweight.getZSliceAt(index)
                u = numpy.where(weight==0,u,u/weight)
            print "u: mean: %.2f std: %.2f" %(numpy.mean(u),numpy.std(u))


            for iz in range( self.nz ):

                Z = self.Zmin + iz

                mat = self.bckprjImageAtZ( u, self.prjmap[index], Z=Z )

                sys.stdout.write( '%i..' %Z )
                sys.stdout.flush()

                volumes[iz] = mat + volumes[iz]

            sys.stdout.write("\n")

        u = numpy.asarray(volumes)
        
        print "Volume: mean: %.2f     std: %.2f    min:  %.2f" %(numpy.mean(u),numpy.std(u),numpy.min(u))
        
        if crossvol==None: u = numpy.abs(numpy.min(u)) + u

        print numpy.min(u)

        for iz in range( self.nz ):
            #u = numpy.asarray(volumes[iz]).astype('float')
            if crossvol!=None: u[iz] = u[iz]*fcross.getZSliceAt(iz)
            mrcVolume.setZSliceAt( iz, u[iz] )
            
        mrcVolume.updateHeader()

        if showVol: os.system("imod %s" %(output))


def createWeightVolume( vol1, vol2, output="weight.mrc"):
    '''Create the weigting volume'''

    v1 = mrc.MRCFile(vol1)
    v2 = mrc.MRCFile(vol2)

    weight = mrc.MRCFile(output)
    weight.setHeader(v1.nx,v1.ny,v1.nz)

    v_1 = []
    v_2 = []

    for iz in range(v1.nz):
        v_1.append(v1.getZSliceAt(iz))
        v_2.append(v2.getZSliceAt(iz))
    
    v_1 = numpy.asarray(v_1).astype('float')
    v_2 = numpy.asarray(v_2).astype('float')

#    alpha = (v_2-numpy.mean(v_2))/4.8/numpy.std(v_2)    # for A
    alpha = (v_2-numpy.mean(v_2))/7.0/numpy.std(v_2)    # for B
#    alpha = numpy.where(alpha>=1.0,1.0,alpha)
#    alpha = numpy.where(alpha<=-1.0,-1.0,alpha)
    weight_ = numpy.exp(-alpha)


    for iz in range(v1.nz):
        weight.setZSliceAt(iz,weight_[iz])

    weight.updateHeader()


def mixVol2( f1_name, f2_name, limits, output_name="test-mix.mrc", med_filt_window=5, method="norm2" ):
    '''Given two reconstruction of the same, generate a mixture of both,
    trying to minimize artifacts'''

    f1 = mrc.MRCFile(f1_name)
    f2 = mrc.MRCFile(f2_name)

    # Get the offsets in x and y

    if limits==None:
        limits = {}
        limits[f1_name] = numpy.array([[0,f1.nx],[0,f1.ny]])
        limits[f2_name] = numpy.array([[0,f2.nx],[0,f2.ny]])

    xmin = max(limits[f1_name][0,0],limits[f2_name][0,0])
    xmax = min(limits[f1_name][0,1],limits[f2_name][0,1])
    ymin = max(limits[f1_name][1,0],limits[f2_name][1,0])
    ymax = min(limits[f1_name][1,1],limits[f2_name][1,1])

    width = xmax - xmin
    height = ymax - ymin

    start_x1 = xmin - limits[f1_name][0,0]
    start_y1 = ymin - limits[f1_name][1,0]
    start_x2 = xmin -limits[f2_name][0,0]
    start_y2 = ymin - limits[f2_name][1,0]

    end_x1 = start_x1 + width
    end_y1 = start_y1 + height
    end_x2 = start_x2 + width
    end_y2 = start_y2 + height

    f = mrc.MRCFile(output_name)
    f.setHeader(width-2,height-2,f1.nz-2)

    # Helper function --- begin---
    def delta( u1, u2, method):

        if method=="norm1":
            delta1 = numpy.abs(u1[0,1:-1,2:]-u1[0,1:-1,1:-1]) + numpy.abs(u1[0,1:-1,:-2]-u1[0,1:-1,1:-1])
            delta1 = delta1 + numpy.abs(u1[0,2:,1:-1]-u1[0,1:-1,1:-1]) + numpy.abs(u1[0,:-2,1:-1]-u1[0,1:-1,1:-1])
            delta1 = delta1 + numpy.abs(u1[0,1:-1,1:-1]-u1[1,1:-1,1:-1]) + numpy.abs(u1[0,1:-1,1:-1]-u1[-1,1:-1,1:-1])
            delta2 = numpy.abs(u2[0,2:,1:-1]-u2[0,1:-1,1:-1]) + numpy.abs(u2[0,:-2,1:-1]-u2[0,1:-1,1:-1])
            delta2 = delta2 + numpy.abs(u2[0,1:-1,2:]-u2[0,1:-1,1:-1]) + numpy.abs(u2[0,1:-1,:-2]-u2[0,1:-1,1:-1])
            delta2 = delta2 + numpy.abs(u2[0,1:-1,1:-1]-u2[1,1:-1,1:-1]) + numpy.abs(u2[0,1:-1,1:-1]-u2[-1,1:-1,1:-1])

        if method=="norm2":
            delta1 = numpy.abs(u1[0,1:-1,2:]-u1[0,1:-1,1:-1])**2 + numpy.abs(u1[0,1:-1,:-2]-u1[0,1:-1,1:-1])**2
            delta1 = delta1 + numpy.abs(u1[0,2:,1:-1]-u1[0,1:-1,1:-1])**2 + numpy.abs(u1[0,:-2,1:-1]-u1[0,1:-1,1:-1])**2
            delta1 = delta1 + numpy.abs(u1[0,1:-1,1:-1]-u1[1,1:-1,1:-1])**2 + numpy.abs(u1[0,1:-1,1:-1]-u1[-1,1:-1,1:-1])**2
            delta2 = numpy.abs(u2[0,2:,1:-1]-u2[0,1:-1,1:-1])**2 + numpy.abs(u2[0,:-2,1:-1]-u2[0,1:-1,1:-1])**2
            delta2 = delta2 + numpy.abs(u2[0,1:-1,2:]-u2[0,1:-1,1:-1])**2 + numpy.abs(u2[0,1:-1,:-2]-u2[0,1:-1,1:-1])**2
            delta2 = delta2 + numpy.abs(u2[0,1:-1,1:-1]-u2[1,1:-1,1:-1])**2 + numpy.abs(u2[0,1:-1,1:-1]-u2[-1,1:-1,1:-1])**2

        return delta1, delta2
    # Helper function --- end ---
    
    for iz in range(1,f1.nz-1):
        print "Tilt #: %i/%i" %(iz,f1.nz-2)
        u1 = (f1.getZSliceAt(iz),f1.getZSliceAt(iz+1),f1.getZSliceAt(iz-1))
        u2 = (f2.getZSliceAt(iz),f2.getZSliceAt(iz+1),f2.getZSliceAt(iz-1))
        u1 = numpy.array(u1)
        u2 = numpy.array(u2)
        u1 = u1[:,start_x1:end_x1,start_y1:end_y1]
        u2 = u2[:,start_x2:end_x2,start_y2:end_y2]
        u1 = (u1-numpy.mean(u1))/numpy.std(u1)
        u2 = (u2-numpy.mean(u2))/numpy.std(u2)
        delta1, delta2 = delta( u1, u2, method)
        u1 = u1[0,1:-1,1:-1]
        u2 = u2[0,1:-1,1:-1]
        x = delta1/(delta1+delta2)
        x = scipy.signal.medfilt2d(x,med_filt_window)
        #x = 0.5                 #######################
        u = (1-x)*u1 + x*u2
        f.setZSliceAt(iz-1,u)

    f.updateHeader()


if __name__ == '__main__':


    def mixtureLikeCrossVal():

        directory="/home/sphan/data/crosstest"

        f1_name = "fhv6a_z_-50.0.out"
        f2_name = "fhv6b_z_-50.0.out"

        f1_name = os.path.join(directory,f1_name)
        f2_name = os.path.join(directory,f2_name)

        limits = { f1_name:numpy.array([[0,398],[0,399]]), f2_name:numpy.array([[-16,403],[-19,402]]) }
        limits = None
        
        mixVol2( f1_name, f2_name, limits )


    def crossReconstruction( directory="/home/sphan/data/crosstest", basename="fhv6" ):

        import txbr

        project = txbr.TxBRproject( directory, basename )
        project.load()

        limits = [ project.reconstruction.origin.z, project.reconstruction.end.z ]

        projection = project.series[0].projection

        prjmap1 = {}

        for index in range(projection.numberOfExposures()):
            prjmap1[index] = projection.getAffineTransformation(index)


        stack1 = txbr.utilities.MRCStack("%s.st" %project.series[0].basename)
        reco1 = CrossReconstruction( stack1 )
        reco1.setRotationAxis( project.series[0].rotAxis )
        reco1.setProjectionMap( prjmap1 )
#        reco1.filter()
#        reco1.bckprj( limits=limits, output=project.series[0].getReconstructionFileName(), showVol=True)
#
#
        projection = project.series[1].projection

        prjmap2 = {}

        for index in range(projection.numberOfExposures()):
            prjmap2[index] = projection.getAffineTransformation(index)

        stack2 = txbr.utilities.MRCStack("%s.st" %project.series[1].basename)
        reco2 = CrossReconstruction( stack2 )
        reco2.setRotationAxis( project.series[1].rotAxis )
        reco2.setProjectionMap( prjmap2 )
#        reco2.filter()
#        reco2.bckprj( limits=limits, output=project.series[1].getReconstructionFileName(), showVol=True)
#

#        createWeightVolume( project.series[0].getReconstructionFileName(), project.series[0].getReconstructionFileName(),output="weightA.mrc")
        createWeightVolume( project.series[1].getReconstructionFileName(), project.series[1].getReconstructionFileName(),output="weightB.mrc")
#
#        txbr.prj.project(project,0,0, srcVolume="weightA.mrc",output="proj0to0.mrc")
        txbr.prj.project(project,1,1, srcVolume="weightB.mrc",output="proj1to1.mrc")
#
#        reco1.bckprj( limits=limits, output="testA.mrc", crossvol="weightA.mrc", weight="proj0to0.mrc", showVol=True)
        reco2.bckprj( limits=limits, output="testB.mrc", crossvol="weightB.mrc", weight="proj1to1.mrc", showVol=True)

        mixVol2( "testA.mrc", "testB.mrc", None, output_name="test.mrc" )


    # Different options

    options = { 1:mixtureLikeCrossVal, 2:crossReconstruction }

    options[1]()
#    options[1]( basename = "FHV102709-29" )
#    options[2]( )