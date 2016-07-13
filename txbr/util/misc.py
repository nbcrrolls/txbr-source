import numpy
import scipy.spatial
import util


def connectPairWise( pts1, pts2, tf, dm=5.0, verbose=False ):
    '''Compare two sets of (3d or 2d) points given a transformation *tf*.
    *tf*(*pts2*) = *pts1*.

    :param pts1: The first set of points
    :param pts2: The second set of points
    :param tf: The transformation that roughly maps the two sets with
         *pts1=tf(pts2)*
    '''

    in2D = False
    in3D = False

    if isinstance(tf,util.PolynomialTransformation3D): in3D = True
    if isinstance(tf,util.PolynomialTransformation2D) or isinstance(tf,util.AffineTransformation2D): in2D = True

    if in2D:
        pts1_ = pts1[:,:2].copy()
        pts2_ = pts2[:,:2].copy()
        pts2_ = tf.forward(pts2_)

    if in3D:
        pts1_ = pts1[:,:3].copy()
        pts2_ = pts2[:,:3].copy()
        pts2_ = tf.forward(pts2_)

        pts1_ = pts1_[:,:2]
        pts2_ = pts2_[:,:2]

    # Try to match points between the two models

    tree1 = scipy.spatial.KDTree(pts1_)
    tree2 = scipy.spatial.KDTree(pts2_)

    d12,i12 = tree1.query(pts2_)
    d21,i21 = tree2.query(pts1_)

    if dm!=None:
        check1 = (i12[i21]==range(len(pts1_)))*(d12[i21]<dm)
        check2 = (i21[i12]==range(len(pts2_)))*(d21[i12]<dm)
    else:
        check1 = (i12[i21]==range(len(pts1_)))
        check2 = (i21[i12]==range(len(pts2_)))

    index1 = numpy.where(check1==True)
    index2 = i21[index1]

    index1 = numpy.squeeze(index1)

    d = d21[index1]

    n = index1.size

    if verbose:
        print "Found %i matching points!" %(n)
        print "Min/Max distance between two matching points: [%.1f,%.1f]" %(numpy.min(d),numpy.max(d))

    return index1,index2


def clusterize( pts, r=5.0, min_cluster_size=1 ):

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

    return numpy.row_stack(res)



if __name__ == "__main__":
    print "Misc module"
