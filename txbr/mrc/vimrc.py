import pylab
import matplotlib
import mrc

def view(filename,x=None,y=None,z=1):

    print z

    file = mrc.MRCFile(filename)

    if x==None: x = file.nx/2
    if y==None: y = file.ny/2
    if z==None: z = file.nz/2


    u_x = file.getXSliceAt(x)
    u_y = file.getYSliceAt(y)
    u_z = file.getZSliceAt(z)

    kws = {}
    if file.mode!=16:
        kws['cmap'] = matplotlib.cm.gray

    pylab.figure()

    pylab.imshow(u_z,**kws)

#    pylab.subplot(221)
#    pylab.imshow(u_x,**kws)
#
#    pylab.subplot(223)
#    pylab.imshow(u_z,**kws)
#
#    pylab.subplot(224)
#    pylab.imshow(u_y,**kws)

    pylab.show()

