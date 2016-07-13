"""
The **quickalign** module has been designed to allow for automatic
alignment in Electron Tomography. It mainly consists of the single class **QuickAlign**
whose purpose is to generate 3D estimates for some marker positions, given a set
of tracks on the 2D electron micrographs and a set of rough projection maps.

The asumption behind the **quickalign** module is:

- to have a good initial estimate for the projection map,O

- to have beads that can be used for alignment.

Projection maps can be calculated first by using *fiducialless* methods for instance
(see :py:mod:txbr.onthefly.quickreco).
"""

import os.path
import multiprocessing
import copy
import numpy
import numpy.linalg
import scipy.spatial
import scipy.spatial.distance
import scipy.cluster.vq
import scipy.stats
import txbr
import modl
import util

from txbr.txbrconf import DATA_DIR
from txbr.txbrconf import DETECT_DIR
from txbr.txbrconf import ALIGN_DIR

from txbr import log

GLOBAL_ALIGN = "global-alignement"
LOCAL_ALIGN = "local-alignement"

STRAIGHTFORWARD_CLUSTER_METHOD = "straighforward-cluster-distance"
HIERARCH_AGG_CLUSTER_METHOD = "hierarchical-agglomerative"
DD_CLUSTER_METHOD = "direct-distance"
K_MEANS_CLUSTER_METHOD = "k-means"

def __proj_markers_on_slice__( prj, XYZ, index, tracks, r ):
    """
    Helper function used in the tracking procedure (for multiprocessing purpose). It
    allows to calculate the 3D marker projections on a given tilt specified by its *index*.

    :param index: Index of the projection slice
    :param tracks: self.points[index]
    """

    pts = prj.forward(XYZ)
    pts[:,2] = index

    if len(tracks)==0: # There is no reference projections

        pts_ = numpy.empty_like(pts)
        pts_[:,:] = numpy.NaN

        i_ = numpy.empty_like(pts[:,0])
        i_[:] = numpy.NaN

        print "Tilt #%i:  0 markers " %(index)

    else:

        eps = numpy.inf
        niter = 50
        iter = 0

        while eps>0.01 and iter<niter: # Do a center of mass iterative correction

            iter += 1

            pts1_ = pts
            pts2_ = tracks

            tree1 = scipy.spatial.KDTree(pts1_[:,:2])
            tree2 = scipy.spatial.KDTree(pts2_[:,:2])

            d12,i12 = tree1.query(pts2_[:,:2])
            d21,i21 = tree2.query(pts1_[:,:2])

            check1 = (i12[i21]==range(len(pts1_)))
            check2 = (i21[i12]==range(len(pts2_)))

            index1 = numpy.where(check1==True)
            index2 = i21[index1]

            diff = pts2_[index2] - pts1_[index1]

            mix = 0.5

            eps = numpy.linalg.norm(numpy.median(diff,axis=0))
            pts[:,:2] = pts[:,:2] + mix*numpy.median(diff,axis=0)[:2]

        npick = numpy.sum(check1)

        # Add the absolute distance constraint

        res_err = numpy.sqrt(numpy.sum(diff**2,axis=1))
        
        mean_ = numpy.mean(res_err)
        std_ = numpy.std(res_err)
        mix_ = 1.0

        thrshold = min(r,mean_+mix_*std_)

        check1 = check1*(d21<=thrshold)
        check2 = check2*(d12<=thrshold)

        index1 = numpy.where(check1==True)
        index2 = i21[index1]

        pts_ = numpy.empty_like(pts)
        pts_[:,:] = numpy.NaN
        pts_[index1] = pts2_[index2]

        i_ = numpy.empty((pts_.shape[0]))
        i_[:] =  numpy.NaN
        i_[index1] = index1

        nkeep = numpy.sum(check1)

        rate = float(nkeep)/npick

        print "Tilt #{0:3}:  picked/kept={1:3}/{2:3} [<{3:.1f} (max:{4:3.1f})]  [eps={5:.2e}] [from pool of {6:}]".format(index,npick,nkeep,thrshold,r,eps,len(tracks))
        if rate<0.5:
            print "\tdiff: med: {0:s}   mean: {1:s}  std: {2:s}".format(str(numpy.median(diff,axis=0)[:2]),str(numpy.mean(diff,axis=0)[:2]),str(numpy.std(diff,axis=0)[:2]))

        return pts_,i_,pts


def __proj_markers_on_slice_star__( args ):
    """
    Wrapper to the helper function __proj_markers_on_slice__ used for multiprocessing purpose.
    """

    return __proj_markers_on_slice__( *args )



class QuickAlign:
    """A class that implements a quick alignment"""

    def __init__( self, basename, points, projmap, xmin, xmax, ymin, ymax, dist2d=None,
                  directory=".", workDirectory=".", bin=None, method=GLOBAL_ALIGN, saveClusterModel=False):

        self.basename = basename
        self.directory = directory
        self.workDirectory = workDirectory
        self.bin = bin
        self.points = points
        self.projmap = projmap
        self.projmap_init = copy.deepcopy(projmap)
        self.saveClusterModel = saveClusterModel

        if self.bin==None: self.bin = 1

        # define the bounding box
        
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

#        wx = 0.5*(self.xmax -self.xmin)
#        wy = 0.5*(self.ymax -self.ymin)
#
#        self.xmin -= wx
#        self.xmax += wx
#        self.ymin -= wy
#        self.ymax += wy

        self.zmin = -max(ymax-ymin,xmax-xmin)/3.0
        self.zmax = -self.zmin

        print "[xmin,xmax]=[%.2f, %.2f]   [ymin,ymax]=[%.2f, %.2f]     [zmin,zmax]=[%.2f, %.2f]" \
                        %(self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)

        #  --- Convergence treshold parameters ---

        # self.nmax is the maximum number of crossed points that will be considered
        # self.crossDistanceThreshold is the maximum tolerated distance between two backprojection
        # lines corresponding to a cross-point.
        # self.distanceThresholdIn2D is the maximum tolerated distance between a projected point
        # estimate and its measure
        # self.minimumClusterSize is the minimum number of points to form a 3D cluster

        if method==LOCAL_ALIGN:
            self.nmax = None
            self.crossDistanceThreshold = 2.5
            self.distanceThresholdIn2D = 2.5
        elif  method==GLOBAL_ALIGN:
            self.nmax = 150
            self.crossDistanceThreshold = 4.0
            self.minimumClusterSize = 3.0
            self.distanceThresholdIn2D = 3.0


        #  ---------------------------------------

        if dist2d!=None:
            self.distanceThresholdIn2D = dist2d

        for index in range(len(self.points)):
            if len(self.points[index])!=0:
                self.points[index] = numpy.row_stack(self.points[index])

        # Redefinition of the projection parameters

        T,M,N = [],[],[]

        for index,p in self.projmap.iteritems():
            T.append(p.T)
            M.append(p.M)
            N.append(numpy.cross(p.M[0,:],p.M[1,:]))

        self.T = numpy.row_stack(T)
        self.M = numpy.array(M)
        self.N = numpy.array(N) # Vector director of the projections


    def saveTracks( self, tracks, lc=0.5, doShortenTracks=True ):
        """Routine to save the tracks as a model file in the alignment directory.

        :param tracks: A list containing the traces as numpy arrays.
        :param lc: The trace length ratio below which a trace is not saved if variable doShortenTracks
            set to true
        :param doShortenTracks: Traces that are to short wont be saved
        """

        n = len(self.projmap)

        if doShortenTracks:
            lengths = [ len(trace) for trace in tracks ]
            ntracks_init = len(tracks)
            tracks = [ trace for trace in tracks if len(trace)>lc*n ]
            ntracks_final = len(tracks)
            log.info("Number of tracks decreased from {0:} to {1:}".format(ntracks_init,ntracks_final))



        if self.bin==None or self.bin==1:  # In case, there is no binning
            mrc_file = os.path.join( self.directory, "%s.st" %self.basename )
        else:
            mrc_file = os.path.join( self.workDirectory, DATA_DIR, "bin%i" %self.bin, "%s.st" %self.basename )

        marker_file = os.path.join( self.workDirectory, ALIGN_DIR, "bin%i" %self.bin, "%s.mrk" %self.basename )

        modl.saveAs3DModel( marker_file, tracks , mrcFileName=mrc_file, doShow=False )


    def backCross( self, pts1, pts2, d1, d2, xmin, xmax, ymin, ymax, zmin, zmax, r=2.0, nmax=200, bounds=[True,True,True] ):
        '''
        d1 and d2: two numpy arrays of shape (2,4) defining the projection maps
        '''

        if nmax==None:
            nmax = max(len(pts1),len(pts2))

        if len(pts1)>=nmax:
            pts1 = pts1[:nmax].copy()
        else:
            pts1 = pts1.copy()

        if len(pts2)>=nmax:
            pts2 = pts2[:nmax].copy()
        else:
            pts2 = pts2.copy()

        a1 = numpy.cross(d1[0,1:],d1[1,1:])    # Vector director of projection 1
        a2 = numpy.cross(d2[0,1:],d2[1,1:])    # Vector director of projection 2

        n1 = len(pts1)
        n2 = len(pts2)

        M = numpy.zeros((6,6))

        M[0,:3] = d1[0,1:]
        M[1,:3] = d1[1,1:]
        M[2,3:] = d2[0,1:]
        M[3,3:] = d2[1,1:]
        M[4,:3] = a1
        M[4,3:] = -a1
        M[5,:3] = a2
        M[5,3:] = -a2

        Minv = numpy.linalg.pinv(M)

        pts1[:,0] -= d1[0,0]
        pts1[:,1] -= d1[1,0]

        pts2[:,0] -= d2[0,0]
        pts2[:,1] -= d2[1,0]

        H1 = numpy.tensordot(Minv[:,0:2],pts1[:,:2],axes=((1),(1)))
        H1 = numpy.repeat(H1,n2,axis=1)
        H1.resize((6,n1,n2))

        H2 = numpy.tensordot(Minv[:,2:4],pts2[:,:2],axes=((1),(1)))
        H2 = numpy.repeat(H2,n1,axis=0)
        H2.resize((6,n1,n2))

        H = H1 + H2

        H1 = H[:3,:,:]
        H2 = H[3:,:,:]
        H = (H1+H2)/2.0

        if H.size==0:
            return numpy.zeros((0,3)),numpy.zeros((0,2))

        D = H2 - H1
        D = numpy.apply_along_axis(numpy.linalg.norm, 0, D)  # distance

        accept = numpy.ones_like(D)
        accept = (D<=r)
        if bounds[0]: accept = accept*(H[0,:,:]>xmin)*(H[0,:,:]<xmax)
        if bounds[1]: accept = accept*(H[1,:,:]>ymin)*(H[1,:,:]<ymax)
        if bounds[2]: accept = accept*(H[2,:,:]>zmin)*(H[2,:,:]<zmax)

        neighborConstraint = True # For a symmetrical neighbor constraint

        if neighborConstraint:

            D = numpy.where(accept,D,numpy.Inf)

            i12 = numpy.argmin(D,axis=0)
            i21 = numpy.argmin(D,axis=1)

            check1 = i12[i21]==range(len(pts1))
            check2 = i21[i12]==range(len(pts2))

            index1 = numpy.where(check1==True)
            index2 = i21[index1]

            index1 = numpy.squeeze(index1)
            index2 = numpy.squeeze(index2)

            points = H[:,index1,index2]

            return points.T,numpy.row_stack((index1,index2)).T

        i1,i2 = numpy.where(accept)
        points = H[:,i1,i2]

        return points.T,numpy.row_stack((i1,i2)).T


    def localAlignement( self, lc=0.5, nc=250, loadLinkFiles=False ):
        """
        This routine implements the tracking by analyzing the detected tracks correspondances
        in a sequential manner. A given track on one view is associated to another in the following
        view, if the corresponding backprojection rays cross-over within a distance less than the
        criteria *self.crossDistanceThreshold*.

        :param lc: The trace length ratio below which a track won't be considered further in
            the alignement.
        :param nc: The maximum number of tracks that will be kept once this local alignement
            routine has been implemented. Used when calling the self.reduce routine.
        :param loadLinkFiles: Load the linked markers from a prevous session if it is available.
        """

        n = len(self.projmap)
        views = range(n)

        step = 1
        distances = []

        for index in views:
            index1 = index
            index2 = index + step
            if index2>=n: break
            linkfile = os.path.join( self.workDirectory, DETECT_DIR, "bin%i" %self.bin, "%s.links-local.%03i.txt" %(self.basename,index))
            if os.path.exists(linkfile) and loadLinkFiles:
                D = numpy.loadtxt(linkfile)
                D = numpy.resize(D,(D.size/2,2))
            else:
                pts1 = self.points[index1]
                pts2 = self.points[index2]
                K1 = self.projmap[index1].M[:2,:]
                K2 = self.projmap[index2].M[:2,:]
                T1 = self.projmap[index1].T[:2]
                T2 = self.projmap[index2].T[:2]
                d1 = numpy.column_stack((T1,K1))
                d2 = numpy.column_stack((T2,K2))
                H,D = self.backCross( pts1, pts2, d1, d2, self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax,
                                    r=self.crossDistanceThreshold, nmax=self.nmax, bounds=[True,True,True] )
                print "Observation from crossing views %3i and %3i: %i points" %(index1,index2,H.shape[0])
                numpy.savetxt(linkfile,D)
            distances.append(D)

        self.connectTracks( distances, lc=lc, nc=nc )


    def connectTracks( self, D, lc=0.5, nc=None, eps_err=50 ):
        """
        Tracks are built from the pairwise relationships between tracks from sequential views.

        :param D: A list containing links between tracks from two sequential views.
        :param lc: The trace length ratio below which a track won't be considered further in
            the alignement.
        :param nc: The maximum number of tracks that will be kept once this local alignement
            routine has been implemented. Used when calling the self.reduce routine.
        """

        doRemoveProblematicClusters = True
        doShortenChains = True
        doReduceChains = True
        doConnectClusters = False

        n = len(self.projmap)

        cluster_id = 0

        cluster_map = {}    # Associates a couple (view,track) to a cluster
        clusters = []

        # Initialize the cluster mapping

        for index,pts in enumerate(self.points):
            for i in range(len(pts)):
                key = "%i-%i" %(index,i)
                cluster_map[key] = []

        # Step (i): create the clusters from the knowledge of sequential links

        for index,lnks in enumerate(D):

            for i,j in lnks:

                key1 = "%i-%i" %(index,i)
                key2 = "%i-%i" %(index+1,j)

                if not cluster_map[key1] and not cluster_map[key2]: # Make a new cluster

                    clusters.append(self.createCluster({ index:i, index+1:j }))
                    cluster_map[key1].append(cluster_id)
                    cluster_map[key2].append(cluster_id)
                    cluster_id += 1

                elif cluster_map[key1] and cluster_map[key2]:

                    # We skip the case where both markers already belong to two clusters
                    continue

                elif cluster_map[key1]:    # Otherwise add to an existing cluster

                    ids1 = cluster_map[key1]

                    if len(cluster_map[key1])>1: # Only pick one cluster which *1* belongs to
                        d = numpy.inf
                        for id in ids1:
                            d_ = clusters[id].distanceFromPoint(self.points[index+1][j])
                            if numpy.max(d_)<d:
                                d = d_[-2]
                                c_id = id
                        ids1 = [c_id]

                    # Store the information

                    cluster_map[key2].extend(ids1)

                    for id in ids1:
                        clusters[id].dict[index+1] = j
    
        print "Distribution of length for all the chains"

        self.lengthDistributionInfo(clusters)

        # Step (ii): Deal with the problematic clusters ?

        if doRemoveProblematicClusters:

            problematic_clusters = []

            for index,cluster in enumerate(clusters):
                if not cluster.validate(eps=eps_err,verbose=True):
                    problematic_clusters.append(index)

            npb = len(problematic_clusters)
            ncluster = len(clusters)

            print "Number of problematic cluster (err>%.1f): %i/%i" %(eps_err,npb,ncluster)

            clusters_ids_ = [  index for index in range(ncluster) if not index in problematic_clusters ]
            clusters_ = [  clusters[index] for index in clusters_ids_ ]

        # Only consider the chains over a certain length lc*n

        if doShortenChains:

            large_chains, large_chains_cluster_id = self.removeSmallClusters( clusters_, clusters_ids_, lc*n )

            nlrgcluster = len(large_chains)
            ncluster_ = len(clusters_)

            print "Number of chain longer than %i: %i/%i" %(lc*n, nlrgcluster, ncluster_)

            self.lengthDistributionInfo(large_chains)

            clusters_ids_ = large_chains_cluster_id
            clusters_ = large_chains

        # Step (iii) Locate the 3D marker location

        XYZ = self.locateClusterPoints( clusters_ )

        # Step (iv): eventually reduce the number of cluster to nc

        if doReduceChains:

            clusters_, clusters_ids_, XYZ, index = self.reduce( clusters_, clusters_ids_, XYZ, nc )

            print "Distribution of length for the picked chains"

            self.lengthDistributionInfo(clusters_)

        # Step (v): assemble some clusters

        if doConnectClusters:   # Might need some tweaking

            finals = []
            finals_id = []

            tracks, cluster_sets = self.track( XYZ, r=self.distanceThresholdIn2D )

            cluster_sets = numpy.array(cluster_sets).T  # (nXYZ,nview)

            for index,cluster4view in enumerate(cluster_sets):

               # print "index %i:  cluster %i" %(index,clusters_ids_[index])

                ref_id = clusters_ids_[index]
                ref = clusters[ref_id]

                u  = []
                for i in range(len(cluster4view)):
                    if not numpy.isnan(cluster4view[i]):
                        u.extend(cluster_map["%i-%i" %(i,int(cluster4view[i]))])
                u = set(u)

                currents = []
                ll = ref.length()   # Size of the chain
                for id in u:
                    if id in problematic_clusters: continue
                    current = clusters[id]
                    currents.append(current)
                    if id in problematic_clusters: print "Pb with %i" %id
                    if ref.isClusterCompatibleWith(current) and not id in problematic_clusters: ll += current.length()
                   # print "%5i-%5i: %6s     [%s] [%s]  d=%.2f  %i" %(ref_id,id,ref.isClusterCompatibleWith(current),ref,current,ref.distanceFromCluster(current),ll)
                if len(currents)==0: continue
                final = ref.agglomerate(currents,eps_dist=15.0)
                
                finals_id.append(ref_id)
                finals.append(final)

                clusters_ids_ = finals_id
                clusters_ = finals

                self.lengthDistributionInfo(clusters_)

        # Step (vi): save the data (traces on the tilt series)

        tracks = [ cluster.points() for cluster in clusters_ ]

        self.saveTracks( tracks, doShortenTracks=False )
 

    def globalAlignement( self, lc=0.25, nc=100, ntotc=2000, dist2d=None ):
        """Global alignment routine

        This routine implements the tracking by analyzing the detected tracks correspondances
        globally via a cluster analysis. The self.nmax backprojected rays of the more predominant
        peaks of a view are compared with the ones of another view (10 degrees apart). The potential
        crossover point then lead to cluster that could potentially be a real marker. Their whole tracks
        are calculated and can be used for alignent.

        :param lc: The trace length ratio below which a track won't be considered further in
            the alignement.
        :param nc: The maximum number of tracks that will be kept once this global alignement
            routine has been implemented. Used when calling the self.reduce routine.
        :param ntotc: The maximum number of backprojected cross-points used to determine the 
            relevant clusters
        :param dist2d: The maximum tolerated distance between a projected track and a possible
            measure of it. Assigned to self.distanceThresholdIn2D

        """

        doShortenTracks = True # Would be equivalent so set lc==0

        if dist2d!=None:
            self.distanceThresholdIn2D = dist2d

        n = len(self.projmap)
        half_n = int(n/2)

        # Step (i): Define the rays use for the backprojected cross-points

        views = []

        def ray_asc():
            step = half_n/10
            return range( 0, half_n, step)

        def ray_desc():
            step = half_n/10
            return range( 0, half_n, step)

        def ray_mix():
            relative_rays = []
            nlevel = 5
            for ilevel in range(2,nlevel):
                for i in range(ilevel):
                    indexOfRay = i*half_n/ilevel
                    if not indexOfRay in relative_rays:
                        relative_rays.append(indexOfRay)
            return relative_rays

        options = { 1:ray_asc, 2:ray_desc, 3:ray_mix }

        relative_rays = options[3]()

        for index in relative_rays:
            if index==0: continue
            views.append((half_n,half_n+index))
            views.append((half_n,half_n-index))

        # Step (ii): determine the possible candidate for the 3D track locatons from this set of rays.

        XYZ = {}    # The dictionary containing the crossing points for each view pair
        links = {}   # The dictionary containing the corresponding marker indices in both view for each crossing points
        ntot = 0

        for index,(index1,index2) in enumerate(views):
            if index1<0 or index1>=n: break
            if index2<0 or index2>=n: break
            key = "%i-%i" %(index1,index2)
            pts1 = self.points[index1]
            pts2 = self.points[index2]
            K1 = self.projmap[index1].M[:2,:]
            K2 = self.projmap[index2].M[:2,:]
            T1 = self.projmap[index1].T[:2]
            T2 = self.projmap[index2].T[:2]
            d1 = numpy.column_stack((T1,K1))
            d2 = numpy.column_stack((T2,K2))
            H,D = self.backCross( pts1, pts2, d1, d2, self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax,
                                r=self.crossDistanceThreshold, nmax=self.nmax )
            print "Observation from crossing views %3i and %3i: %i points (between [%i,%i])" %(index1,index2,H.shape[0],len(pts1),len(pts2))
            XYZ[key] = H
            links[key] = D
            linkfile = os.path.join( self.workDirectory, DETECT_DIR, "bin%i" %self.bin, "%s.links-global.%03i-%03i.txt" %(self.basename,index1,index2))
            numpy.savetxt(linkfile,D)
            ntot += len(H)
            if index%2==1 and ntot>ntotc: break # We dont want too many points to clusterize

        print "Total number of cross-observations: %s" %(len(XYZ))

        # Step (iii): Identify the clusters. Reduce their numbers if needed

        #min_cluster_size = min(0.5*len(XYZ),self.minimumClusterSize)
        min_cluster_size = min(0.35*len(XYZ),self.minimumClusterSize)
        clusters, cluster_ids, _3dpts = self.clusterize( XYZ, links, min_cluster_size=min_cluster_size )

        clusters, cluster_ids, _3dpts, index = self.reduce( clusters, cluster_ids, _3dpts, nc )

        # Step (iv): Identify the traces in the tilt series and eventually keep the longer ones

        tracks, cluster_sets = self.track( _3dpts, r=self.distanceThresholdIn2D )

        # to do: populate the clusters
            
#        if doShortenTracks:
#            tracks = [ trace for trace in tracks if len(trace)>lc*n ]
#
#        if self.saveClusterModel:
#            modl.saveAs3DModel( "tata.mod", _3dpts, int(xmax), int(ymax), int(zmax), doShow=False )
#
#        marker_file = os.path.join( self.workDirectory, ALIGN_DIR, "bin%i" %self.bin, "%s.mrk" %self.basename )
#        mrc_file = os.path.join( self.workDirectory, DATA_DIR, "bin%i" %self.bin, "%s.st" %self.basename )
#
#        modl.saveAs3DModel( marker_file, tracks , mrcFileName=mrc_file, doShow=False )
#
        self.saveTracks( tracks, lc=lc, doShortenTracks=doShortenTracks )

        return tracks


    def clusterize( self, XYZ, links, method=HIERARCH_AGG_CLUSTER_METHOD, r=8.0, min_cluster_size=4, nmax=None ):
        '''Consider a XY mesh.

        :paramXYZ: A list containing numpy arrays of points
        :param links:
        :param method: The method use for this clusterization process. It can be
            - STRAIGHTFORWARD_CLUSTER_METHOD
            - HIERARCH_AGG_CLUSTER_METHOD
            - DD_CLUSTER_METHOD
            - K_MEANS_CLUSTER_METHOD
        :param r: Represents the maximum size of of cluster
        :param min_cluster_size: The minimum number of points in a cluster
        :param nmax: The maximum number of points in XYZ that would be use for processing the clusters.
            Could be None or an integer. If nmax is a finite number, it will the constraint will be spread
            out evenly between each pair of views.
        :returns: A 3-tuple with the list of clusters, the list of cluster ids and a numpy array with the
            locations of the clusters
        '''

        pts = []
        if nmax!=None:  # Reduce the number of points for performance issue
            nmean = nmax/len(XYZ)
            for key,points in XYZ.iteritems():
                pts.append(points[:nmean])
        else:
            for key,points in XYZ.iteritems():
                pts.append(points[:])

        pts = numpy.squeeze(numpy.row_stack(pts))

        if method==STRAIGHTFORWARD_CLUSTER_METHOD:

            pts = {}
            tiltref = None  # We work with a reference tilt

            for key,link in links.iteritems():
                tilt = key.split("-")
                tilt1 = int(tilt[0])   # All the tilt1 should be the same
                tilt2 = int(tilt[1])
                if tiltref==None:
                    tiltref = tilt1
                elif tilt1!=tiltref:
                    continue
                for index,l in enumerate(link):
                    if not l[0] in pts:
                        pts[l[0]]= [XYZ[key][index]]
                    else:
                        pts[l[0]].append(XYZ[key][index])

            for index,points in pts.iteritems():
                std = numpy.std(points,axis=0)
                print "Index #%i  npoints=%i  std=%s" %(index,len(points),std)
                print points

            import sys
            sys.exit(0)
                

        elif method==HIERARCH_AGG_CLUSTER_METHOD:

            distances = scipy.spatial.distance.pdist(pts)
            Y = scipy.cluster.hierarchy.linkage(distances)
            T = scipy.cluster.hierarchy.fcluster(Y, r, criterion='distance')

            clusters = numpy.unique(T)
            ncluster = len(clusters)

            dict = {}
            res = []
            cst = min_cluster_size

            while len(res)==0 and cst!=0:

                print "Loop: len(res)=%i    cst=%i" %(len(res),cst)

                for icluster in range(ncluster):
                    _pts_ = pts[numpy.argwhere(T==clusters[icluster])]
                    if len(_pts_)<cst:
                        dict[clusters[icluster]] = None
                    else:
                        dict[clusters[icluster]] = numpy.mean(_pts_,axis=0)
                        res.append(numpy.mean(_pts_,axis=0))

                cst -= 1

            if cst==0:
                log.warning("No real cluster aggregate have been found.")

            res = numpy.row_stack(res)

            ncluster = len(res)

            cluster_ids = range(ncluster)
            clusters = [ self.createCluster({}) for id in cluster_ids ]
            
            print "A total of %i clusters have been found." %(ncluster)

            return cluster_ids, clusters, res

        elif method==DD_CLUSTER_METHOD:

            tree = scipy.spatial.KDTree( pts )
            nearest = tree.query_ball_point( pts, r )

            length = numpy.vectorize(len)
            acc_pts = length(nearest)

            sorted_index = numpy.argsort(acc_pts)

            itemindex = numpy.where(acc_pts[sorted_index]==min_cluster_size)
            itemindex = numpy.min(itemindex)

            sorted_index = sorted_index[itemindex:]

            final = []
            redundant = []

            for index in sorted_index:
                if not index in redundant:
                    final.append(index)
                redundant.extend(nearest[index])

            print "A total of %i clusters have been found." %(len(final))

            return pts[final]

        elif method==K_MEANS_CLUSTER_METHOD:

            ncluster = 100

            codebook, distortion = scipy.cluster.vq.kmeans(pts[sorted_index],ncluster,thresh=D)

            return codebook

        else:

            return numpy.zeros((0,3))

    def __map_regression__( self, projmap, order=3, projmap0=None ):
        """Performs a regression on the set of projection maps, or its residual from
        an initial set *projmap0*"""

        tilts = []
        T,M = [],[]

        for index,p in projmap.iteritems():
            tilts.append(index)
            if projmap0:
                p0 = projmap0[index]
                T.append(p.T-p0.T)
                M.append(p.M-p0.M)
            else:
                T.append(p.T)
                M.append(p.M)

        tilts = numpy.asarray(tilts)
        T = numpy.asarray(T)
        M = numpy.asarray(M)

        shape0 = M.shape
        M.resize((shape0[0],M.size/shape0[0]))

        T_reg = numpy.empty_like(T)
        M_reg = numpy.empty_like(M)

        for index in range(T.shape[1]):
            pfit = numpy.polyfit(tilts,T,order)
            T_reg[:,index] = numpy.poly1d(pfit[:,index])(tilts)

        for index in range(M.shape[1]):
            pfit = numpy.polyfit(tilts,M,order)
            M_reg[:,index] = numpy.poly1d(pfit[:,index])(tilts)

        M_reg.resize(shape0)

        for index,p in projmap.iteritems():
            
            if projmap0:
                p0 = projmap0[index]
                K = numpy.column_stack((T_reg[index]+p0.T,M_reg[index]+p0.M))
            else:
                K = numpy.column_stack((T_reg[index],M_reg[index]))
                
            projmap[index]= util.AffineTransformation3D(K)

        return projmap


    def track( self, XYZ, r=3.0, pickFromDetectedMarkers=True, adjustProjectionMap=False, adjust3DTracks=False ):
        """Track the detected beads on the micrograph series. The procedure consists
        in projecting the 3D position of the beads throughout the entire series. The
        closest detected mark less than a distance r is then considered the corresponding
        projection. The final traces are stored in a file with extension ".mrk".

        :param XYZ: A numpy array containing the cluster model, that is the hypothetical
            position of the bead markers.
        :param pickFromDetectedMarkers: If True (default value) only pick the tracks from the
            detected markers. Otherwise refine take the track from the inferred positons.
        :param r: The maximum tolerance distance to assign a (detected) bead to a track
        :param adjustProjectionMap: Refine the projection map from tracks and 3D marker positions
            when more than 4 objects are available on a tilt micrograph.
        :returns: A 2-tuple containing the list of point-track positions and a list containing
        their corresponding view indices as a numpy array. For the latter one, if a track was
        not found for a 3D point at a given tilt, the index is assigned to numpy.NaN.
        """

        # The loop only kicks in if the 3D marker locations are adjusted

        adjustProjectionMapLocally = True
        
        self.projmap = self.__map_regression__(self.projmap,projmap0=self.projmap_init) # smooth out the projection maps globally

        XYZ = numpy.asarray(XYZ)

        ntilt = len(self.projmap)
        ntrack = XYZ.shape[0]
        tracks = []

        cluster_sets = []

        if not pickFromDetectedMarkers:

            for index,prj in self.projmap.iteritems():

                pts = prj.forward(XYZ)
                pts[:,2] = index

                tracks.append( pts )

        else:

            loop = 0
            ntot0 = -1
            ntot = 0

            while ntot>ntot0:

                loop += 1

                tracks = []
                residual_error = []
                cluster_sets = []

                pool = multiprocessing.Pool(processes=txbr.PROCESS_POOL_SIZE)

                tasks = [[ prj, XYZ, index, self.points[index] , r ] for index,prj in self.projmap.iteritems() ]
                results = pool.imap(__proj_markers_on_slice_star__,tasks)

                for res in results:
                    pts_,i_,pts = res
                    tracks.append( pts_ )
                    residual_error.append(pts_-pts)
                    cluster_sets.append(i_)

                pool.close()
                pool.terminate()

                # Adjust the 3D track positions

                if adjust3DTracks:

                    T = numpy.zeros((2*ntilt))
                    M = numpy.zeros((2*ntilt,3))

                    for index,prj in self.projmap.iteritems():
                        T[2*index] = prj.T[0]
                        M[2*index,:] = prj.M[0,:]
                        T[2*index+1] = prj.T[1]
                        M[2*index+1,:] = prj.M[1,:]

                    tracks_ = numpy.dstack(tracks)
                    tracks_ = numpy.rollaxis( tracks_, 2, 1)

                    XYZnew = numpy.zeros_like(XYZ)

                    for index,trace in enumerate(tracks_):
                        trace_ = trace[:,:2].flatten()
                        ii = numpy.where(~numpy.isnan(trace_))
                        trace_ = trace_[ii] - T[ii]
                        Minv_ = numpy.linalg.pinv(M[ii])
                        XYZnew[index] = numpy.squeeze(numpy.dot(Minv_,trace_))
                        print "index# %i (%s)  %s  ->  %s " %(index,len(trace_)/2,XYZ[index],XYZnew[index])

                    XYZ = XYZnew

                # Refine the projection map globally by shifting the rotation center

                if adjustProjectionMap:

                    T = numpy.zeros((2*ntilt))
                    M = numpy.zeros((2*ntilt,3))

                    for index,prj in self.projmap.iteritems():
                        T[2*index] = prj.T[0]
                        M[2*index,:] = prj.M[0,:]
                        T[2*index+1] = prj.T[1]
                        M[2*index+1,:] = prj.M[1,:]

                    residual_error_ = numpy.dstack(residual_error)
                    residual_error_ = numpy.rollaxis( residual_error_, 2, 1)
                    residual_error_ = residual_error_[:,:,:2].flatten()

                    alpha = numpy.repeat(M,ntrack,axis=0)

                    ii = numpy.where(~numpy.isnan(residual_error_))

                    residual_error_ = residual_error_[ii]
                    alpha_inv_ = numpy.linalg.pinv(alpha[ii])

                    C_ = numpy.dot(alpha_inv_,residual_error_)

                    print "Axis mvt: %s" %(str(C_))

                    mix = 5.0

                    for index,prj in self.projmap.iteritems():

                        prj.T[0] += mix*(C_[0]-numpy.dot(prj.M[0,:],C_))
                        prj.T[1] += mix*(C_[1]-numpy.dot(prj.M[1,:],C_))


                if adjustProjectionMapLocally: # Refine individually the projection map for each view

                    for index,prj in self.projmap.iteritems():

                        B = tracks[index][:,:2].flatten()

                        M = numpy.zeros((2*ntrack,8))
                        M[::2,0] = 1
                        M[::2,1] = XYZ[:,0]
                        M[::2,2] = XYZ[:,1]
                        M[::2,3] = XYZ[:,2]
                        M[1::2,4] = 1
                        M[1::2,5] = XYZ[:,0]
                        M[1::2,6] = XYZ[:,1]
                        M[1::2,7] = XYZ[:,2]

                        ii = numpy.where(~numpy.isnan(B))

                        if len(ii)<4: continue

                        P = numpy.dot(numpy.linalg.pinv(M[ii]),B[ii])

                        prj.T[0] = P[0]
                        prj.T[1] = P[4]
                        prj.M[0,:] = P[1:4]
                        prj.M[1,:] = P[5:8]

                if not adjust3DTracks:
                    break
                else:
                    ntot0 = ntot
                    ntot = 0.0
                    for x in tracks:
                        ntot += numpy.sum(~numpy.isnan(x))/3.0/len(XYZ)
                    print "Loop: %i  Track length=%.1f" %(loop,ntot)


        # Finalize the traces by reshuffling the marker list.
        # Eliminate indices where no beads are found

        tracks = numpy.dstack(tracks)
        tracks = numpy.rollaxis( tracks, 2, 1)

        tracks = list(tracks)
    
        for index,t in enumerate(tracks):   # for each marker...
            iii = numpy.where(~numpy.any(numpy.isnan(t),axis=1))
            tracks[index] = t[iii]

        self.lengthDistributionInfo(tracks)
      
        return tracks,cluster_sets


    def locateClusterPoints ( self, clusters ):
        """Locate 3D cluster points from a set of clusters within the reconstruction

            :param clusters: A list of clusters
        """

        pts = [ cluster.locate3DClusterPoint() for cluster in clusters if cluster.length()>0 ]

        if len(pts)>0:
            return numpy.row_stack([ cluster.locate3DClusterPoint() for cluster in clusters if cluster.length()>0 ])
        else:
            return numpy.zeros((0,3))



    def removeSmallClusters( self, clusters, cluster_ids, lc=20, verbose=True ):
        """
        Remove the cluster shorter than lc

        :param clusters: A list of clusters
        :param cluster_ids: The corresponding list of cluster ids
        :param XYZ: A numpy array of shape (n,3)
        :param ncluster: The maximum number of points that are returned.
        """

        large_chains = []
        large_chains_cluster_id = []

        for cluster,id in zip(clusters,cluster_ids):
            if cluster.length()>lc:
                large_chains.append(cluster)
                large_chains_cluster_id.append(id)
                if verbose: print "Keep cluster #%i: %s" %(id,cluster)

        return large_chains,large_chains_cluster_id


    def reduce( self, clusters, cluster_ids, XYZ, ncluster=300, verbose=False ):
        """
        Reduce a set of points to only *ncluster* points if needed. This routine
        tries to keep a good distribution of points by eliminating the ones that
        are parts of big clusters.

        :param clusters: A list of clusters
        :param cluster_ids: The corresponding list of cluster ids
        :param XYZ: A numpy array of shape (n,3)
        :param ncluster: The maximum number of points that are returned.
        :returns: A 4-tuple containing a ndarray (shape (n.1)) of the final clusters,
            a ndarray (shape (n.1)) of the final cluster ids, a ndarray (shape (n.3)) of
            the final cluster position and the list of index corresponding to the final
            clusters. The parameter *n* is the final number of clusters (:math:`n \leq n_{cluster}`)
        """

        if ncluster==None:
            return clusters, cluster_ids, XYZ, range(len(XYZ))

        clusters = numpy.asarray(clusters)
        cluster_ids = numpy.asarray(cluster_ids)

        if len(clusters)<ncluster:
            return clusters, cluster_ids, XYZ, range(len(clusters))

        XYZ = numpy.asarray(XYZ)

        codebook,distortion = scipy.cluster.vq.kmeans(XYZ,ncluster)

        ncluster = codebook.shape[0]    # The returned number of cluster might be slightly different

        code,dist = scipy.cluster.vq.vq(XYZ, codebook)

        index = -1*numpy.ones((ncluster),dtype='int')

        for i,c in enumerate(code):
            if index[c]<0 or dist[i]<dist[index[c]]:
                index[c] = i
                if verbose: print "Keep cluster #%i: %s" %(cluster_ids[i],clusters[i])

        clusters_red_ = clusters[index]
        cluster_ids_red_ = cluster_ids[index]

        print "Number of clusters was reduced from %i to %i" %(len(clusters),len(clusters_red_))

        return clusters_red_,cluster_ids_red_,XYZ[index],index


    def lengthDistributionInfo( self, tracks ):
        """Some statistics on the length distribution
        
        :param tracks: A list containing numpy arrays of shape (n,3) corresponding to a
            the projection trace for each marker.
        """

        ntilt = len(self.points)
        
        L = []
        for trace in tracks:
            if type(trace)==numpy.ndarray:
                L.append(len(trace))
            if isinstance(trace,self.__Cluster):
                L.append(trace.length())

        print "Total number of cluster: %i" %(len(tracks))

        nbins = 10
        hist,bin_edges = numpy.histogram(L,nbins,range=(0.5,ntilt+0.5))

        for index,h in enumerate(hist):
            print "%i chains of size in [%.1f,%.1f]" %( h, bin_edges[index], bin_edges[index+1] )

        Lmin = numpy.min(L)
        Lmax = numpy.max(L)
        Lmean = numpy.mean(L)

        print "Mean Track length: %f for %i patches  [min=%i,max=%i]" %(Lmean,len(tracks),Lmin,Lmax)

        return L


    def createCluster( self, dict ):

        return QuickAlign.__Cluster( self, dict )


    class __Cluster(object):
        """An inner class for easily taking into account traces"""

        def __init__( self, qa, dict ):

            self.qa = qa    # The QuickAlign object
            self.dict = dict


        def points( self ):
            """Marker points that belong to this cluster"""

            pts = []

            for index,track in self.dict.iteritems():
                pts.append(self.qa.points[index][track])

            if len(pts)==0:
                return numpy.zeros((0,3))

            return numpy.row_stack(pts)


        def length( self ):
            """The number of markers within this cluster"""

            return len(self.dict)


        def min( self ):
            """The number of markers within this cluster"""

            return min(self.dict.keys())


        def max( self ):
            """The number of markers within this cluster"""

            return max(self.dict.keys())


        def locate3DClusterPoint( self ):
            """Calculate the position of the cluster point within the reconstruction"""

            XYZ,Err = self.estimateClusterPoint()

            return XYZ


        def estimateClusterPoint( self ):
            """Calculate the position of the cluster point within the reconstruction"""

            n = len(self.dict)

            index = numpy.zeros((n),dtype='int')
            x = numpy.zeros((n),dtype='int')
            y = numpy.zeros((n),dtype='int')

            i = 0
            for view,track in self.dict.iteritems():
                x[i],y[i],index[i] = self.qa.points[view][track]
                i += 1

            index = index.astype('int')

            B = numpy.zeros((2*n))
            B[:n] = x[:] - self.qa.T[index,0]
            B[n:] = y[:] - self.qa.T[index,1]

            M_ = numpy.zeros((2*n,3))
            M_[:n,0] = self.qa.M[index,0,0]
            M_[:n,1] = self.qa.M[index,0,1]
            M_[:n,2] = self.qa.M[index,0,2]
            M_[n:,0] = self.qa.M[index,1,0]
            M_[n:,1] = self.qa.M[index,1,1]
            M_[n:,2] = self.qa.M[index,1,2]

            Minv = numpy.linalg.pinv(M_)

            XYZ = numpy.dot(Minv,B)

            Err =  numpy.dot(M_,XYZ)-B
            Err = numpy.sqrt(Err[:n]**2 + Err[n:]**2)

            return XYZ,Err


        def validate( self, eps=20.0, verbose=False ):
            """Calculate the position of the cluster point within the reconstruction"""

            XYZ,Err = self.estimateClusterPoint()

            if verbose: print "Cluster 3D point ambiguity: Err %.1f for %i points" %(numpy.sum(Err),len(self.dict))

            return numpy.max(Err)<eps


        def removeOutliers(self):

            XYZ,Err = self.estimateClusterPoint()

            mean = numpy.mean(Err)
            std = numpy.std(Err)

            indexOfOutliers = numpy.where(Err>5.0*std)

            print "outliers: %i" %(len(self.dict))
            print Err
            print "Mean: %.1f  std: %.1f" %(numpy.mean(Err),numpy.std(Err))
            print indexOfOutliers


        
        def isProblematic( self ):

            s = set(self.dict.keys())

            return len(s)!=len(self.dict.keys())


        def isClusterCompatibleWith( self, cluster ):
            """Check is this cluster is compatible with another one"""

            s1 = set(self.dict.keys())
            s2 = set(cluster.dict.keys())
            
            return s1.isdisjoint(s2)


        def distanceFromCluster( self, cluster ):

            X1 = self.locate3DClusterPoint()
            X2 = cluster.locate3DClusterPoint()

            return numpy.linalg.norm(X1-X2)


        def distanceFromPoint( self, point ):

            i0 = int(point[2])
            n0 = self.qa.N[i0]

            M0 = numpy.array([ [ self.qa.M[i0,0,1], self.qa.M[i0,0,2] ] , [ self.qa.M[i0,1,1], self.qa.M[i0,1,2] ] ])
            P0 = numpy.dot( numpy.linalg.pinv(M0), numpy.array([ point[0]-self.qa.T[i0,0], point[1]-self.qa.T[i0,1] ]) )

            distances = []

            for index,track in self.dict.iteritems():

                pt = self.qa.points[index][track]

                i = int(pt[2])
                n = self.qa.N[index]

                v = numpy.cross(n0,n)                
                v = v/numpy.linalg.norm(v)

                M = numpy.array([ [ self.qa.M[i,0,1], self.qa.M[i,0,2] ] , [ self.qa.M[i,1,1], self.qa.M[i,1,2] ] ])
                P = numpy.dot( numpy.linalg.pinv(M0), numpy.array([ pt[0]-self.qa.T[i,0], pt[1]-self.qa.T[i,1] ]) )

                H = numpy.concatenate(([0.0],(P-P0)))

                d = numpy.linalg.norm(numpy.dot(H,v))

                distances.append(d)

            return distances


        def agglomerate( self, clusters, eps_dist=15.0 ):
            """Agglomerate the two clusters"""

            d = []

            for cluster in clusters:
                d.append( self.distanceFromCluster(cluster) )

            indices = numpy.argsort(d)

            agg_clust = self.qa.createCluster( {} )

            for index in indices:
                if not agg_clust.isClusterCompatibleWith(clusters[index]): continue
                if d[index]>eps_dist: continue
                agg_clust.dict.update(clusters[index].dict)

            return agg_clust


        def __str__( self ):

            list_str = ""
            
            for token in sorted(self.dict.keys()):
                if token-1 in self.dict.keys() and token+1 in self.dict.keys():
                    continue
                if not token-1 in self.dict.keys() and token+1 in self.dict.keys():
                    list_str += "%i-" %(token)
                if not token+1 in self.dict.keys():
                    list_str += "%i " %(token)


            return "Cluster: (l=%i)   %s" %( len(self.dict), list_str )

        

if __name__ == '__main__':

    def prj_test_points( directory, basename, prjmap ):

        model_file = os.path.join( directory, "txbr-align", "%s.mod" %basename)

        m = modl.Model(model_file)
        m.loadFromFile()
        
        fid = m.points()

        fid[:,0] += - offset.x
        fid[:,1] += - offset.y
        fid[:,2] += - offset.z

        points = []

        for key,map in prjmap.iteritems():
            pj = map.forward(fid)
            pj[:,2] = key
            points.append(pj)

        return points


    option1 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "model":"FHV-18a.bead.mod", "basename":"FHV-18a"  }
    option2 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "model":"titi.fid", "basename":"FHV-18a"  }
    option3 = { "directory":"/ncmirdata3/sphan/otf/mon", "model":"SCN_night_OTO_cellchain001a.bead.mod", "basename":"SCN_night_OTO_cellchain001a"  }
    option4 = { "directory":"/ncmirdata3/sphan/otf/mon2", "model":"SCN_night_OTO_cellchain001a.bead.mod", "basename":"SCN_night_OTO_cellchain001a"  }
    option5 = { "directory":"/ncmirdata3/sphan/otf/25", "model":"FHV-25-2a.bead.mod", "basename":"FHV-25-2a"  }

    options = option5

    directory = options["directory"]
    basename = options["basename"]
    fmodel = options["model"]

    # Load the detected beads.

    m = modl.Model(fmodel)
    m.loadFromFile()

    obj =  m.objects[0]

    points = []

    for cont in obj.contours:
        points.append(cont.points)

    xmin = m.xoffset
    ymin = m.yoffset
    xmax = m.xmax
    ymax = m.ymax

    points = modl.loadPoints(fmodel)

    # Open the TxBR project to load the projection maps.

    project = txbr.TxBRproject( directory, basename )
    project.load()
    project.series[0].shiftPrealignTransform(eps=-1)  # Shift to st based pjmaps

    rotation = project.reconstruction.getRotation()
    offset = project.reconstruction.getOffset()

    origin = project.reconstruction.origin
    end = project.reconstruction.end

    xmin,ymin,zmin = origin.x,origin.y,origin.z
    xmax,ymax,zmax = end.x,end.y,end.z

    print "x:%f:%f  x:%f:%f  x:%f:%f" %(xmin,xmax,ymin,ymax,zmin,zmax)

    projection = project.series[0].projection
    
    prjmap = {}

    for index in range(projection.numberOfExposures()):
        prjmap[index] = projection.getAffineTransformation(index)



#    points = prj_test_points( directory, basename, prjmap )

    # Perform the quick alignment

    qa = QuickAlign( basename, points, prjmap, xmin, xmax, ymin, ymax )
    
