import numpy,pylab
import modl

def test( input, base, nx=None, ny=None, title=None ):

    X,Y = numpy.array_split(input,2,axis=1)
    U,V = numpy.array_split(base-input,2,axis=1)

    print title + ':'
    print '  Deviation Range in X: [%f,%f]' %(numpy.amin(U),numpy.amax(U))
    if nx!=None:
        print '                        [%.2f%%,%.2f%%]' %(numpy.amin(U)*100.0/nx,numpy.amax(U)*100.0/nx)
    print '  Deviation Range in Y: [%f,%f]' %(numpy.amin(V),numpy.amax(V))
    if ny!=None:
        print '                        [%.2f%%,%.2f%%]' %(numpy.amin(U)*100.0/ny,numpy.amax(U)*100.0/ny)

    pylab.figure()
    pylab.quiver(X,Y,U,V)
    pylab.title(title)


if __name__ == '__main__':

    #model = modl.Model('/Users/sph/Electron_Tomography/txbr/data/NaokoStack/Naoko3A_stack_refine.mod')
    model = modl.Model('/Users/sph/Naoko3A_stack_refine.mod')
    model.loadFromFile()

    nx = model.xmax
    ny = model.ymax

    for index,interface in enumerate(model.objects):

        pts = [contour.points for contour in interface.contours if contour.points.shape[0]==2]
        pts = numpy.row_stack(pts)
        input = pts[::2,:2]
        base = pts[1::2,:2]

        test(input,base,nx=nx,ny=ny,title='Interface %i' %index)

    pylab.show()

