import numpy
import mrc

def undecimate(file,factor=2):

    f = mrc.MRCFile(file)
    nx = f.nx
    ny = f.ny
    nz = f.nz

    fout = mrc.MRCFile(file+'undec',template=file)
    fout.setHeader(factor*nx,factor*ny,factor*nz,mode=f.mode)

    v = numpy.zeros((factor*nx,factor*ny,3))

    for iz in range(nz):
        u = f.getZSliceAt(iz)
        for ichannel in range(3):
            for ix in range(factor):
                for iy in range(factor):
                    xp = 1.0*ix/factor
                    yp = 1.0*iy/factor
                    print 'ix=%i  iy=%i  xp=%f  yp=%f' %(ix,iy,xp,yp)
                    v[ix:-factor:factor,iy:-factor:factor,ichannel] = (1.0-xp)*(1.0-yp)*u[:-1,:-1,ichannel] + \
                                                        xp*(1.0-yp)*u[1:,:-1,ichannel] + \
                                                        (1.0-xp)*yp*u[:-1,1:,ichannel] + \
                                                        xp*yp*u[1:,1:,ichannel]

        fout.setZSliceAt(iz*factor,v)

    for iz in range(nz-1):
        u1 = fout.getZSliceAt(iz*factor)
        u2 = fout.getZSliceAt((iz+1)*factor)
        for i in range(1,factor):
            zp = 1.0*i/factor
            v= (1.0-zp)*u1 + zp*u2
            fout.setZSliceAt(iz*factor+i,v)

    u1 = fout.getZSliceAt((nz-1)*factor-1)
    u2 = fout.getZSliceAt((nz-1)*factor)
    for i in range(1,factor):
        zp = float(i)
        v = (1.0+zp)*u2 - zp*u1
        fout.setZSliceAt(nz*factor-factor+i,v)


if __name__ == '__main__':
    undecimate('/Users/sph/Electron_Tomography/txbr/dataset/LM/Slice2/Slice2.bin16.mrc',factor=2)
