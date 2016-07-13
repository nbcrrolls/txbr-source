import os
import os.path
import sys
import re
import numpy
import pylab
import scipy.integrate
import mrc

def crossValidate( directory=".", basenames=(), mode = 'outlier-max', h=1.0, doPlot=False ):
    '''
    Variable mode: in [ 'int-simps', 'int-romb', 'int-regular', 'outlier-max', 'outlier-min', 'most-probable', 'dequantization', 'remove' ]
    Variable h: used when variable mode is set to "dequantization"
    '''

    verbose = False

    print "directory: %s" %directory
    print "basenames: %s" %basenames

    if not os.path.exists(directory): return
    if len(basenames)==0: return

    if mode=='remove': h = int(h)

    file_out = { 'int-simps':"output.int.simps.mrc",
                 'int-romb':"output.int.romb.mrc",
                 'int-regular':"output.int.reg.mrc",
                 'outlier-max':"output.outlier.max.mrc",
                 'outlier-min':"output.outlier.min.mrc",
                 'most-probable':"output.most.probable.mrc",
                 'variance':"output.variance.mrc",
                 'std-normalize':"output.std.norm.mrc",
                 'dequantization':"output.dequant.%.2f.mrc" %h,
                 'remove':"output.remove.%i.mrc" %h }

    file_out = os.path.join( directory, file_out[mode] )

    outliers_flag = 0.0   # ( >0: remove positive outliers, <0: remove negative outliers, 0: both)
    outlier_threshold = 1.5 # ( if 0.0: nothing removed )

    outlier_threshold = numpy.abs(outlier_threshold)    # Make sure it is positive

    files = os.listdir(directory)

    tilts = {}
    input = []

    for basename in basenames:
        tilts[basename] = []
        for file in files:
            match = re.match(basename + '*.(\d*)$', file)
            if not match: continue
            print "%s: %i" %(basename,int(match.group(1)))
            tilts[basename].append(int(match.group(1)))
        tilts[basename].sort()
        input.extend([ os.path.join(directory,"%s.%i" %(basename,tilt)) for tilt in tilts[basename] if os.path.exists(os.path.join(directory,"%s.%i" %(basename,tilt)))])

    ntilt = len(input)

    if (ntilt==0): sys.exit()

    nxstart,nystart,nzstart = [],[],[]
    nxstop,nystop,nzstop = [],[],[]

    for fin in input:
        f = mrc.MRCFile(fin)
        nxstart.append( f.nxstart )
        nystart.append( f.nystart )
        nzstart.append( f.nzstart )
        nxstop.append( f.nxstart + f.nx )
        nystop.append( f.nystart + f.ny )
        nzstop.append( f.nzstart + f.nz )

    nxstart = numpy.round(numpy.max(nxstart)).astype(int)
    nystart = numpy.round(numpy.max(nystart)).astype(int)
    nzstart = numpy.round(numpy.max(nzstart)).astype(int)
    nxstop = numpy.round(numpy.min(nxstop)).astype(int)
    nystop = numpy.round(numpy.min(nystop)).astype(int)
    nzstop = numpy.round(numpy.min(nzstop)).astype(int)

    nx = nxstop - nxstart
    ny = nystop - nystart
    nz = nzstop - nzstart

    print "(nxstart,nystart,nzstart)=(%i,%i,%i)   (nx,ny,nz)=(%i,%i,%i)" %(nxstart,nystart,nzstart,nx,ny,nz)

    if nx<=0 or ny<=0 or nz<=0:return

    for fin in input: print fin

    crop = [5,2]
    crop = numpy.min(numpy.vstack((crop,[len(tilts)/2,len(tilts)/2])),axis=0)

    z = range(nzstart,nzstop)
 #   z = [12,50]

    indices = [(222,223,224),(243,256,256)]

    fout = mrc.MRCFile(file_out)
    fout.setHeader(nx,ny,len(z))

    # Do some "cross validation" between tilts/series

    for index,iz in enumerate(z):
        
        u = []
        
        for file_in in input:
            print "Open file %s at index %i" % (file_in, index)
            fin = mrc.MRCFile(file_in)
            z1 = numpy.round(fin.nzstart).astype(int)
            z2 = z1 + fin.nz
            if z1>iz or iz>=z2:
                continue
            slice = fin.getZSliceAt(iz-z1)
            x1 = nxstart - numpy.round(fin.nxstart).astype(int)
            x2 = x1 + nx
            y1 = nystart - numpy.round(fin.nystart).astype(int)
            y2 = y1 + ny
            slice = slice[x1:x2,y1:y2]
            u.append(slice)

        u = numpy.array(u)
        u = numpy.sort(u, axis=0)
        m = numpy.mean(u, axis=0)
        std = numpy.std(u, axis=0)

        if outliers_flag>0: outliers = (u-m)/std > outlier_threshold
        elif outliers_flag<0: outliers = (u-m)/std < -outlier_threshold
        else: outliers = numpy.abs(u-m)/std > outlier_threshold

        if verbose: print "Outliers: %.3f%%" %(float(numpy.sum(outliers))/outliers.size)

        if crop[1]!=0: u = u[crop[0]:-crop[1],:,:]
        elif crop[0]!=0: u = u[crop[0]:,:,:]

        t = numpy.rollaxis(u,0,3)[indices]

        #u = numpy.log(u - numpy.min(u) + 0.1)

        if mode=='int-simps': u = scipy.integrate.simps(u,axis=0)
        elif mode=='int-romb': u = scipy.integrate.romb(u,axis=0)
        elif mode=='outlier-max': u = numpy.max(u,axis=0)
        elif mode=='outlier-min': u = numpy.min(u,axis=0)
        elif mode=='remove':
            v = numpy.sort(numpy.abs(u-m),axis=0)
            index_ = v.shape[0] - h
            if index_<0: index_=0
            print "%i/%i" %(index_,v.shape[0])
            threshold = v[index_,:,:]
            u = numpy.where(v<=threshold,u,0.0)
            #u = numpy.where(v<=threshold,0.0,u)
            u = numpy.sum(u,axis=0)
        elif mode=='std-normalize':
          #  u = m + (u-m)/std
            u = u/std
            u = numpy.sum(u,axis=0)
        elif mode=='variance':
            u = numpy.std(u, axis=0)
        elif mode=='dequantization':
            u = h*numpy.log(numpy.sum(numpy.exp(u/h),axis=0))
        elif mode=='most-probable':
            print u.shape
            ntilt,nx,ny = u.shape
            w = numpy.zeros((nx,ny))
            for ix in range(nx):
                    for iy in range(ny):
                            der = numpy.zeros(ntilt)
                            der[1:-1] = u[2:,ix,iy]-u[:-2,ix,iy]
                            der[0] = u[1,ix,iy]-u[0,ix,iy]
                            der[-1] = u[-1,ix,iy]-u[-2,ix,iy]
                            inx = numpy.argmin(der)
                            w[ix,iy] = u[inx,ix,iy]
            u = w
        else: u = numpy.sum(u,axis=0)

        # print u.shape

        #u = numpy.log(u - numpy.min(u) + 0.1)
        #u = numpy.sum(u,axis=0) - numpy.min(u) -numpy.max(u)

        fout.setZSliceAt(index,u)


    fout.updateHeader()
    os.system("imod %s" % (file_out))

    if doPlot:
        pylab.figure()
        pylab.plot(t.T)
        pylab.show()

if __name__ == '__main__':

    directory = "/home/sphan/data/artifact-removal/cross-validation/filtered"
   # directory = "/home/sphan/data/artifact-removal/cross-validation/unfiltered"

    basenames = [ "fhv6a_z_-50.0.out", "fhv6b_z_-50.0.out" ]
 #   basenames = [ "fhv6a_z_-50.0.out" ]

##    crossValidate( directory=directory, basenames=basenames, mode = 'outlier-max' )
#    crossValidate( directory=directory, basenames=basenames, mode = 'remove', h=10 )
#    crossValidate( directory=directory, basenames=basenames, mode = 'remove', h=20 )
    #crossValidate( directory=directory, basenames=basenames, mode = 'variance' )
    crossValidate( directory=directory, basenames=basenames, mode = 'std-normalize' )



#    for h in [ -1000.0,-100.0,-10.0,-1.0,1.0,10.0,100.0,1000.0]:
#        crossValidate( directory=directory, basenames=basenames, mode = 'dequantization', h=h )
   