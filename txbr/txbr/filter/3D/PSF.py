from math import pi,sin,cos,atan,floor,ceil,sqrt,exp
from time import time
import numpy
import txbr

def get_TxBR_PSF_r(txbrProject):

    tic = time()

    nx = txbrProject.reconstruction.nx
    ny = txbrProject.reconstruction.ny
    nz = txbrProject.reconstruction.nz

    print nx
    print ny
    print nz

    k_inc = 0.1;
    kmax = sqrt(nx*nx+ny*ny+nz*nz)

    PSF = numpy.zeros((nx,ny,nz),dtype=float)

    trajectory = txbrProject.series[0].trajectory
    numberOfExposures = trajectory.getNumberOfExposures()

    for i in range(numberOfExposures):
        u1 = numpy.array([-1,0,trajectory.X_coefficients[i][3]])
        u2 = numpy.array([0,-1,trajectory.Y_coefficients[i][3]])
        u = numpy.cross(u1,u2)
        u = u/sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
        print i
        for r in numpy.arange(-nz+1,nz-1,k_inc):
            r1 = r*u[0]
            r2 = r*u[1]
            r3 = r*u[2]
            #if r1<0 or r1>=nx-1 or r2<0 or r2>ny-1 or r3<0 or r3>=nz-1: continue
            if abs(r1)>=nx-1 or abs(r2)>ny-1 or abs(r3)>=nz-1: continue
            i1 = floor(r1)
            i2 = floor(r2)
            i3 = floor(r3)
            r1_diff = r1-i1
            r2_diff = r2-i2
            r3_diff = r3-i3
            PSF[i1,i2,i3] += (1-r1_diff)*(1-r2_diff)*(1-r3_diff)
            if i1<nx-1: PSF[i1+1,i2,i3] += r1_diff*(1-r2_diff)*(1-r3_diff)
            if i2<ny-1: PSF[i1,i2+1,i3] +=  (1-r1_diff)*r2_diff*(1-r3_diff)
            if i1<nx-1 and i2<ny-1: PSF[i1+1,i2+1,i3] += r1_diff*r2_diff*(1-r3_diff)
            if i3<nz-1: PSF[i1,i2,i3+1] += (1-r1_diff)*(1-r2_diff)*r3_diff
            if i1<nx-1 and i3<nz-1: PSF[i1+1,i2,i3+1] += r1_diff*(1-r2_diff)*r3_diff
            if i2<ny-1 and i3<nz-1: PSF[i1,i2+1,i3+1] += (1-r1_diff)*r2_diff*r3_diff
            if i1<nx-1 and i2<ny-1 and i3<nz-1: PSF[i1+1,i2+1,i3+1] += r1_diff*r2_diff*r3_diff

#    import pylab
#    pylab.figimage(PSF[:,1,:])
#    pylab.show()

    toc = time()
    print toc-tic,' has elapsed'

#    PSF = PSF/PSF[0,0,0]
#
#    print PSF.min()
#    print PSF.max()
#
#    PSF = PSF + 0.01

    return PSF

def get_TxBR_PSF_k(txbrProject):

    tic = time.time()

    nx = txbrProject.reconstruction.nx
    ny = txbrProject.reconstruction.ny
    nz = txbrProject.reconstruction.nz

    k_inc = 1;
    kmax = sqrt(nx*nx+ny*ny+nz*nz)

    PSF = numpy.zeros((nx,ny,nz),dtype=float)

    trajectory = txbrProject.series[0].trajectory
    numberOfExposures = trajectory.getNumberOfExposures()

    for i in range(numberOfExposures):
        u1 = numpy.array([-1,0,trajectory.X_coefficients[i][3]])
        u2 = numpy.array([0,-1,trajectory.Y_coefficients[i][3]])
        u1 = u1/numpy.linalg.norm(u1)
        u2 = u2/numpy.linalg.norm(u2)
        u2 = u2-u1*u2
        #u3 = numpy.cross(u1,u2)
        print i
        for k1 in numpy.arange(-kmax,kmax,k_inc):
            for k2 in numpy.arange(-kmax,kmax,k_inc):
                k = 1
                k = k1*u1 + k2*u2
#                k1_ = k1*u1[0] + k2*u2[0]
#                k2_ = k1*u1[1] + k2*u2[1]
#                k3_ = k1*u1[2] + k2*u2[2]
                if abs(k[0])>=nx or abs(k[1])>=ny or abs(k[2])>=nz: continue
                i1 = floor(k[0])
                i2 = floor(k[1])
                i3 = floor(k[2])
                k1_diff = k[0]-i1
                k2_diff = k[1]-i2
                k3_diff = k[2]-i3
                PSF[i1,i2,i3] += (1-k1_diff)*(1-k2_diff)*(1-k3_diff)
                if i1<nx-1: PSF[i1+1,i2,i3] += k1_diff*(1-k2_diff)*(1-k3_diff)
                if i2<ny-1: PSF[i1,i2+1,i3] +=  (1-k1_diff)*k2_diff*(1-k3_diff)
                if i1<nx-1 and i2<ny-1: PSF[i1+1,i2+1,i3] += k1_diff*k2_diff*(1-k3_diff)
                if i3<nz-1: PSF[i1,i2,i3+1] += (1-k1_diff)*(1-k2_diff)*k3_diff
                if i1<nx-1 and i3<nz-1: PSF[i1+1,i2,i3+1] += k1_diff*(1-k2_diff)*k3_diff
                if i2<ny-1 and i3<nz-1: PSF[i1,i2+1,i3+1] += (1-k1_diff)*k2_diff*k3_diff
                if i1<nx-1 and i2<ny-1 and i3<nz-1: PSF[i1+1,i2+1,i3+1] += k1_diff*k2_diff*k3_diff

#    import pylab
#    pylab.figimage(PSF[:,1,:])
#    pylab.show()

    toc = time.time()
    print toc-tic,' has elapsed'

    PSF = PSF/PSF[0,0,0]

    print PSF.min()
    print PSF.max()

    PSF = PSF + 0.01

    return PSF

def get_PSFa_2D(nx,ny,angle1,inc,angle2):

    PSF = numpy.zeros((nx,ny),dtype=float)
    rmax = sqrt(nx*nx+ny*ny)
    print rmax

    for angle in numpy.arange(angle1,angle2+1,inc):
    #for angle in numpy.arange(-90,91,inc):

        print angle

        if angle==angle1 or angle==angle2:
            t = 0.5
        else:
            t = 1

        r_inc = 0.1

        for r in numpy.arange(-rmax,rmax,r_inc):
            theta = pi/180*angle
            x = r*sin(theta)
            y = r*cos(theta)
            #if abs(x)>nx/2 or abs(y)>ny/2:
            #if abs(x)>nx-1 or abs(y)>ny-1:
            if abs(x)>100 or abs(y)>100:
                continue
#            if r!=0 and (angle<angle1 or angle>angle2):
#                t = abs(r)
#                #t = exp(-1.0/t*t)/t
#                t = 1/t
#                t = 1
#            else:
#                t = 1
            i = floor(x)
            j = floor(y)
            x_diff = x-i
            y_diff = y-j
            PSF[i,j] += t*(1-x_diff)*(1-y_diff)
            if i<nx-1:
                PSF[i+1,j] += t*x_diff*(1-y_diff)
            if j<ny-1:
                PSF[i,j+1] += t*(1-x_diff)*y_diff
            if i<nx-1 and j<ny-1:
                PSF[i+1,j+1] += t*x_diff*y_diff

    PSF = PSF/PSF[0,0]

    return PSF

def get_PSFa_k_2D(nx,ny,angle1,inc,angle2):

    PSF = get_PSFa_2D(ny,nx,angle1,inc,angle2)

    PSF = PSF + 0.01

    return PSF.transpose()


def get_PSFb_2D(nx,ny,angle1,angle2):

    PSF = numpy.zeros((nx,ny),dtype=float)

    PSF_0 = 0.005

    for i in range(0,nx):
        for j in range(0,ny):
            if i<nx/2:
                x = float(i);
            else:
                x = float(i-nx);
            if j<ny/2:
                y = float(j);
            else:
                y = float(j-ny);
            r = sqrt(x*x+y*y)
            if x>=0 and y==0:
                angle = 90
            elif x<0 and y==0:
                angle = -90
            else:
                angle = 180/pi*atan(x/y)
            if r!=0 and angle>=angle1 and angle<=angle2:
                PSF[i,j] = 1/r + PSF_0
            else:
                PSF[i,j] =  PSF_0

    PSF[0,0] = 1.5

    return PSF

def get_PSFb_k_2D(nx,ny,angle1,angle2):

    PSF = numpy.zeros((nx,ny),dtype=float)

    rmax = sqrt(nx*nx+ny*ny)

    for i in range(0,nx):
        for j in range(0,ny):
            if i<nx/2:
                x = float(i);
            else:
                x = float(i-nx);
            if j<ny/2:
                y = float(j);
            else:
                y = float(j-ny);
            r = sqrt(x*x+y*y)
            if x>=0 and y==0:
                angle = 90
            elif x<0 and y==0:
                angle = -90
            else:
                angle = 180/pi*atan(x/y)
            PSF_0 = atan(r/20)*2/pi
            if r!=0 and r<rmax and (angle<=angle1 or angle>=angle2):
                PSF[i,j] = PSF_0/r
            else:
                PSF[i,j] =  0.01
#            if j==10:
#                print '%i %f  %f  %f' %(i,x,y,angle)

    return PSF