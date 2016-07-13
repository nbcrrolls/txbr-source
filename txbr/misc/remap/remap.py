#!/usr/bin/python

import logging, os.path, time, sys
import numpy, numpy.linalg, pylab, scipy.optimize
import mrc, mrc.mrcfile, modl
import misc.remap.remapSupp

import numpy.random
import util.tester

from misc.remap.remapSupp import NO_INTERPOLATION, BILINEAR_INTERPOLATION, BICUBIC_INTERPOLATION

use_ncg = False
ITERATION_NUMBER = 100000
NPRINT = 1000
iter = 0


numpy.set_printoptions(precision=4)

LOG_CONF = os.path.dirname(__file__) + '/log.conf'
logging.config.fileConfig(LOG_CONF)

log = logging.getLogger()

# Transformation types

TRANSLATIONS = 0
AFFINE_TRANSFORMATIONS = 1
PROJECTIVE_TRANSFORMATIONS_1 = 2
PROJECTIVE_TRANSFORMATIONS_2 = 3

TIFF = 'TIFF'
MRC = 'MRC'
WIDTH = 10

'''
Upper Left  :	quadrant #1
Upper Right :	quadrant #2
Lower Right :	quadrant #3
Lower Left  :	quadrant #4
'''


def prepareStichingModel(directory, basename, ntilt):
	'''This routine creates the MRC file that will be used as a stichng model when
	merging the micrographs of the four quadrants of the NCMIR JEOL microscope 4000#2.
	'''

	# Eventually create MRC directory

	mrc_dir = os.path.join(directory, MRC)

	if not os.path.exists(mrc_dir): os.mkdir(mrc_dir)

	# tranform tif into mrc

	files = []

	for i in range(1,ntilt+1):
		for j in range(1,5):
			file_tif = os.path.join(directory, TIFF, '%s_%03i_cam%i.tif' %(basename,i,j))
			file_mrc = os.path.join(directory, MRC, '%s_%03i_cam%i.tif' %(basename,i,j))
			files.append(file_mrc)
			command = 'tif2mrc %s %s' %(file_tif,file_mrc)
			log.info(command)
			os.system(command)

	nx,ny,nz = mrc.mrcfile.size(file_mrc)

	# Merge the quadrants tiff file onto one single MRC file

	final_mrc = mrc.MRCFile('%s.mrc' %basename)
	final_mrc.setHeader(2*nx+WIDTH,2*ny+WIDTH,ntilt)

	section = numpy.zeros((2*nx+WIDTH,2*ny+WIDTH))

	for i in range(1,ntilt+1):
		log.info('Tilt #%i' %i)
		filename1 = os.path.join(directory, MRC, '%s_%03i_cam1.mrc' %(basename,i))
		f = mrc.MRCFile(filename1)
		section[:nx,ny+WIDTH:] = f.getZSliceAt(0)
		filename2 = os.path.join(directory, MRC,'%s_%03i_cam2.mrc' %(basename,i))
		f = mrc.MRCFile(filename2)
		section[nx+WIDTH:,ny+WIDTH:] = f.getZSliceAt(0)
		filename3 = os.path.join(directory, MRC, '%s_%03i_cam3.mrc' %(basename,i))
		f = mrc.MRCFile(filename3)
		section[nx+WIDTH:,:ny] = f.getZSliceAt(0)
		filename4 = os.path.join(directory, MRC, '%s_%03i_cam4.mrc' %(basename,i))
		f = mrc.MRCFile(filename4)
		section[:nx,:ny] = f.getZSliceAt(0)
		final_mrc.setZSliceAt(i-1,section)

	final_mrc.updateHeader()

	os.system('imod %s' %final_mrc.filename)

	return final_mrc


def quadrant(x,y,nx,ny,width):
	'''Returns the quadrant number.'''

	if x<nx and y<ny: return 4
	elif x>=nx+width and y<ny: return 3
	elif x<nx and y>=ny+width: return 1
	elif x>=nx+width and y>=ny+width: return 2

	return -1


def readPoints(filename):
	'''Read the data points from the stiching model'''

	log.info('Load stitching model %s' %filename)

	f = modl.Model(filename)
	f.loadFromFile()

	xmax = f.xmax
	ymax = f.ymax

	nx = (xmax-WIDTH)/2
	ny = (ymax-WIDTH)/2

	pts = numpy.zeros((0,6))

	for object in f.objects:
		for contour in object.contours:
			if contour.points.shape[0]!=2: continue
			log.debug('Contour: %i' %contour.indexOfContour)
			x1,y1 = contour.points[0,:2]
			q1 = quadrant(x1,y1,nx,ny,WIDTH)-1
			x2,y2 = contour.points[1,:2]
			q2 = quadrant(x2,y2,nx,ny,WIDTH)-1
			log.debug('  x=%f  y=%f  -> quadrant %i' %(x1,y1,q1))
			log.debug('  x=%f  y=%f  -> quadrant %i' %(x2,y2,q2))
			pts = numpy.row_stack((pts,[x1,y1,q1,x2,y2,q2]))

	return (pts,nx,ny)


def translate(pts,delta_x,delta_y):

	pts[:,0] += delta_x
	pts[:,3] += delta_x
	pts[:,1] += delta_y
	pts[:,4] += delta_y


def plot_stitching_offset(pts):

	X = pts[:,0]
	Y = pts[:,1]
	U = pts[:,0] - pts[:,3]
	V = pts[:,1] - pts[:,4]

	pylab.figure()
	pylab.quiver(X,Y,U,V,scale=500)
	pylab.show()


def getPureTranslations(pts):

	M = numpy.zeros((0,8))
	T = numpy.zeros((0,2))

	for point in pts:
		x1,y1,q1 = point[0],point[1],int(point[2])
		x2,y2,q2 = point[3],point[4],int(point[5])
		m = [0]*8
		m[2*q1] = 1
		m[2*q1+1] = 1
		m[2*q2] = -1
		m[2*q2+1] = -1
		t = [0]*2
		t[0] = x2 - x1
		t[1] = y2 - y1
		M = numpy.row_stack((M,m))
		T = numpy.row_stack((T,t))

	M = numpy.row_stack((M,[1]*8))
	T = numpy.row_stack((T,[0]*2))

	Mx = M[:,::2]
	Tx = T[:,0]

	Mxinv = numpy.linalg.pinv(Mx)
	Tx =  numpy.dot(Mxinv,Tx)

	My = M[:,1::2]
	Ty = T[:,1]

	Myinv = numpy.linalg.pinv(My)
	Ty = numpy.dot(Myinv,Ty)

	new_pts = pts.copy()

	q1 = pts[:,2].astype('int')
	q2 = pts[:,5].astype('int')

	new_pts[:,0] += Tx[q1]
	new_pts[:,1] += Ty[q1]

	new_pts[:,3] += Tx[q2]
	new_pts[:,4] += Ty[q2]

	resx = new_pts[:,0] - new_pts[:,3]
	resy = new_pts[:,1] - new_pts[:,4]

	res = numpy.sqrt(resx**2+resy**2)

	res = numpy.array([numpy.min(res),numpy.max(res),numpy.mean(res)])

	Tx_ = numpy.zeros((12))
	Ty_ = numpy.zeros((12))

	Tx_[::3] = Tx
	Ty_[::3] = Ty

	Tx_[1::3] = 1.0
	Ty_[2::3] = 1.0

	return (Tx_,Ty_,new_pts,res)


def getAffineTransformations(pts,nx,ny):

	M = numpy.zeros((0,24))
	T = numpy.zeros((0,2))

	for point in pts:

		x1,y1,q1 = point[0],point[1],int(point[2])
		x2,y2,q2 = point[3],point[4],int(point[5])

		m = [0]*24

		m[6*q1] = 1
		m[6*q1+1] = 1
		m[6*q1+2] = x1	# should be 1
		m[6*q1+3] = x1	# should be 0
		m[6*q1+4] = y1	# should be 0
		m[6*q1+5] = y1	# should be 1

		m[6*q2] = -1
		m[6*q2+1] = -1
		m[6*q2+2] = -x2	# should be 1
		m[6*q2+3] = -x2	# should be 0
		m[6*q2+4] = -y2	# should be 0
		m[6*q2+5] = -y2	# should be 1

		t = [0]*2
		t[0] = 0.0
		t[1] = 0.0

		M = numpy.row_stack((M,m))
		T = numpy.row_stack((T,t))

	M = numpy.row_stack((M,[1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0]))
	T = numpy.row_stack((T,[0]*2))

	# Additional Constraint

	M = numpy.row_stack((M,[0,0,nx,0,0,ny,0,0,nx,0,0,ny,0,0,nx,0,0,ny,0,0,nx,0,0,ny]))
	T = numpy.row_stack((T,[4*nx,4*ny]))

	M = numpy.row_stack((M,[0,0,0,nx,ny,0,0,0,0,nx,ny,0,0,0,0,nx,ny,0,0,0,0,nx,ny,0]))
	T = numpy.row_stack((T,[0,0]))

	# Solve the system

	Mx = M[:,::2]
	Tx = T[:,0]

	Mxinv = numpy.linalg.pinv(Mx)
	Tx =  numpy.dot(Mxinv,Tx)

	My = M[:,1::2]
	Ty = T[:,1]

	Myinv = numpy.linalg.pinv(My)
	Ty = numpy.dot(Myinv,Ty)

	new_pts = pts.copy()

	q1 = pts[:,2].astype('int')
	q2 = pts[:,5].astype('int')

	new_pts[:,0] = Tx[3*q1] + Tx[3*q1+1]*pts[:,0] + Tx[3*q1+2]*pts[:,1]
	new_pts[:,1] = Ty[3*q1] + Ty[3*q1+1]*pts[:,0] + Ty[3*q1+2]*pts[:,1]

	new_pts[:,3] = Tx[3*q2] + Tx[3*q2+1]*pts[:,3] + Tx[3*q2+2]*pts[:,4]
	new_pts[:,4] = Ty[3*q2] + Ty[3*q2+1]*pts[:,3] + Ty[3*q2+2]*pts[:,4]

	resx = new_pts[:,0] - new_pts[:,3]
	resy = new_pts[:,1] - new_pts[:,4]


	res = numpy.sqrt(resx**2 + resy**2)

	res = numpy.array([numpy.min(res),numpy.max(res),numpy.mean(res)])

	return (Tx,Ty,new_pts,res)


def getProjectiveTransformations_1(pts, nx, ny, qfix=[0,1]):
	'''The Remap is performed with a projective transformation, with an isotropic
	scaling change.
	'''

	def residual(T,q_fix,ext):

		Tx = T[0:12]
		Ty = T[12:24]
		lambda_ = T[24:]

		x1,y1,q1 = pts[:,0],pts[:,1],pts[:,2].astype('int')
		x2,y2,q2 = pts[:,3],pts[:,4],pts[:,5].astype('int')

		new_pts = numpy.zeros_like(pts)

		X1 = Tx[3*q1] + Tx[3*q1+1]*x1 + Tx[3*q1+2]*y1
		Y1 = Ty[3*q1] + Ty[3*q1+1]*x1 + Ty[3*q1+2]*y1

		X2 = Tx[3*q2] + Tx[3*q2+1]*x2 + Tx[3*q2+2]*y2
		Y2 = Ty[3*q2] + Ty[3*q2+1]*x2 + Ty[3*q2+2]*y2

		lam_1 = 1.0 + lambda_[2*q1]*x1 + lambda_[2*q1+1]*y1
		lam_2 = 1.0 + lambda_[2*q2]*x2 + lambda_[2*q2+1]*y2

		new_pts[:,0], new_pts[:,1] = X1/lam_1, Y1/lam_1
		new_pts[:,3], new_pts[:,4] = X2/lam_2, Y2/lam_2

		resx = new_pts[:,0] - new_pts[:,3]
		resy = new_pts[:,1] - new_pts[:,4]

		res2 = resx**2 + resy**2
		res = numpy.sqrt(res2)

		res2 = numpy.sum(res2)

		global iter

		line_output = ' #%-8i	Mean Res.: %-8.3f	Max. Res.: %-8.3f	Full Res. Sqr.: %-8.3f' %(iter,numpy.mean(res), numpy.max(res), numpy.sum(res2))

		if iter%NPRINT==0:
			print "\r%s" %line_output
		iter += 1

		for iq in q_fix:

			X = Tx[3*iq] + Tx[3*iq+1]*ext[iq,0] + Tx[3*iq+2]*ext[iq,1]
			Y = Ty[3*iq] + Ty[3*iq+1]*ext[iq,0] + Ty[3*iq+2]*ext[iq,1]
			lam = 1.0 + lambda_[2*iq]*ext[iq,0] + lambda_[2*iq+1]*ext[iq,1]
			lx = X/lam - ext[iq,0]
			ly = Y/lam - ext[iq,1]

			if iq%2==0: res2 += lx**2
			if iq%2==1: res2 += ly**2

		res_ = numpy.array([numpy.min(res), numpy.max(res), res2, numpy.mean(res)])

		return (new_pts,res_,resx,resy)


	def der_residual(T,q_fix,ext):

		Tx = T[0:12]
		Ty = T[12:24]
		lambda_ = T[24:]

		x1,y1,q1 = pts[:,0],pts[:,1],pts[:,2].astype('int')
		x2,y2,q2 = pts[:,3],pts[:,4],pts[:,5].astype('int')

		new_pts = numpy.zeros_like(pts)

		X1 = Tx[3*q1] + Tx[3*q1+1]*x1 + Tx[3*q1+2]*y1
		Y1 = Ty[3*q1] + Ty[3*q1+1]*x1 + Ty[3*q1+2]*y1

		X2 = Tx[3*q2] + Tx[3*q2+1]*x2 + Tx[3*q2+2]*y2
		Y2 = Ty[3*q2] + Ty[3*q2+1]*x2 + Ty[3*q2+2]*y2

		lam_1 = 1.0 + lambda_[2*q1]*x1 + lambda_[2*q1+1]*y1
		lam_2 = 1.0 + lambda_[2*q2]*x2 + lambda_[2*q2+1]*y2

		new_pts[:,0] = X1/lam_1
		new_pts[:,1] = Y1/lam_1

		new_pts[:,3] = X2/lam_2
		new_pts[:,4] = Y2/lam_2

		resx = new_pts[:,0] - new_pts[:,3]
		resy = new_pts[:,1] - new_pts[:,4]

		der_res = numpy.zeros(32)

		indices1,indices2 = [],[]

		for iq in range(4):

			indices1.append(numpy.argwhere(q1==iq))
			indices2.append(numpy.argwhere(q2==iq))

			der_res[3*iq] = 2*numpy.sum(resx[indices1[iq]]/lam_1[indices1[iq]]) \
						   - 2*numpy.sum(resx[indices2[iq]]/lam_2[indices2[iq]])
			der_res[3*iq+1] = 2*numpy.sum(x1[indices1[iq]]*resx[indices1[iq]]/lam_1[indices1[iq]]) \
							 - 2*numpy.sum(x2[indices2[iq]]*resx[indices2[iq]]/lam_2[indices2[iq]])
			der_res[3*iq+2] = 2*numpy.sum(y1[indices1[iq]]*resx[indices1[iq]]/lam_1[indices1[iq]]) \
							 - 2*numpy.sum(y2[indices2[iq]]*resx[indices2[iq]]/lam_2[indices2[iq]])

			der_res[12+3*iq] = 2*numpy.sum(resy[indices1[iq]]/lam_1[indices1[iq]]) \
						   - 2*numpy.sum(resy[indices2[iq]]/lam_2[indices2[iq]])
			der_res[12+3*iq+1] = 2*numpy.sum(x1[indices1[iq]]*resy[indices1[iq]]/lam_1[indices1[iq]]) \
							 - 2*numpy.sum(x2[indices2[iq]]*resy[indices2[iq]]/lam_2[indices2[iq]])
			der_res[12+3*iq+2] = 2*numpy.sum(y1[indices1[iq]]*resy[indices1[iq]]/lam_1[indices1[iq]]) \
							 - 2*numpy.sum(y2[indices2[iq]]*resy[indices2[iq]]/lam_2[indices2[iq]])

			der_res[24+2*iq] = - 2*numpy.sum(resx[indices1[iq]]*new_pts[indices1[iq],0]*x1[indices1[iq]]/lam_1[indices1[iq]]) \
						   	   + 2*numpy.sum(resx[indices2[iq]]*new_pts[indices2[iq],3]*x2[indices2[iq]]/lam_2[indices2[iq]]) \
						       - 2*numpy.sum(resy[indices1[iq]]*new_pts[indices1[iq],1]*x1[indices1[iq]]/lam_1[indices1[iq]]) \
							   + 2*numpy.sum(resy[indices2[iq]]*new_pts[indices2[iq],4]*x2[indices2[iq]]/lam_2[indices2[iq]])

			der_res[24+2*iq+1] = - 2*numpy.sum(resx[indices1[iq]]*new_pts[indices1[iq],0]*y1[indices1[iq]]/lam_1[indices1[iq]]) \
						   	     + 2*numpy.sum(resx[indices2[iq]]*new_pts[indices2[iq],3]*y2[indices2[iq]]/lam_2[indices2[iq]]) \
						         - 2*numpy.sum(resy[indices1[iq]]*new_pts[indices1[iq],1]*y1[indices1[iq]]/lam_1[indices1[iq]]) \
						         + 2*numpy.sum(resy[indices2[iq]]*new_pts[indices2[iq],4]*y2[indices2[iq]]/lam_2[indices2[iq]])

		for iq in q_fix:

			if iq%2==0:
				X = Tx[3*iq] + Tx[3*iq+1]*ext[iq,0] + Tx[3*iq+2]*ext[iq,1]
				lam = 1.0 + lambda_[2*iq]*ext[iq,0] + lambda_[2*iq+1]*ext[iq,1]
				lx = X/lam - ext[iq,0]
				der_res[3*iq] += 2.0*lx/lam
				der_res[3*iq+1] += 2.0*lx*ext[iq,0]/lam
				der_res[3*iq+2] += 2.0*lx*ext[iq,1]/lam
				der_res[24+2*iq] -= 2.0*lx*X*ext[iq,0]/lam**2
				der_res[24+2*iq+1] -= 2.0*lx*X*ext[iq,1]/lam**2

			if iq%2==1:
				Y = Ty[3*iq] + Ty[3*iq+1]*ext[iq,0] + Ty[3*iq+2]*ext[iq,1]
				lam = 1.0 + lambda_[2*iq]*ext[iq,0] + lambda_[2*iq+1]*ext[iq,1]
				ly = Y/lam - ext[iq,1]
				der_res[12+3*iq] += 2.0*ly/lam
				der_res[12+3*iq+1] += 2.0*ly*ext[iq,0]/lam
				der_res[12+3*iq+2] += 2.0*ly*ext[iq,1]/lam
				der_res[24+2*iq] -= 2.0*ly*Y*ext[iq,0]/lam**2
				der_res[24+2*iq+1] -= 2.0*ly*Y*ext[iq,1]/lam**2

		return der_res


	def hess_res(T,q_fix,ext):

		Tx = T[0:12]
		Ty = T[12:24]
		lambda_ = T[24:]

		x1,y1,q1 = pts[:,0],pts[:,1],pts[:,2].astype('int')
		x2,y2,q2 = pts[:,3],pts[:,4],pts[:,5].astype('int')

		new_pts = numpy.zeros_like(pts)

		lam_1 = 1.0 + lambda_[2*q1]*x1 + lambda_[2*q1+1]*y1
		lam_2 = 1.0 + lambda_[2*q2]*x2 + lambda_[2*q2+1]*y2

		new_pts[:,0] = (Tx[3*q1] + Tx[3*q1+1]*x1 + Tx[3*q1+2]*y1)/lam_1
		new_pts[:,1] = (Ty[3*q1] + Ty[3*q1+1]*x1 + Ty[3*q1+2]*y1)/lam_1

		new_pts[:,3] = (Tx[3*q2] + Tx[3*q2+1]*x2 + Tx[3*q2+2]*y2)/lam_2
		new_pts[:,4] = (Ty[3*q2] + Ty[3*q2+1]*x2 + Ty[3*q2+2]*y2)/lam_2

		resx = new_pts[:,0] - new_pts[:,3]
		resy = new_pts[:,1] - new_pts[:,4]

		hess_res = numpy.zeros((32,32))

		for iq in range(4):

			indices1 = numpy.argwhere(q1==iq)
			indices2 = numpy.argwhere(q2==iq)

			hess_res[3*iq,3*iq] = 2.0*numpy.sum(1.0/lam_1[indices1]**2) \
								+ 2.0*numpy.sum(1.0/lam_2[indices2]**2)
			hess_res[3*iq,3*iq+1] = 2.0*numpy.sum(x1[indices1]/lam_1[indices1]**2) \
								  + 2.0*numpy.sum(x2[indices2]/lam_2[indices2]**2)
			hess_res[3*iq,3*iq+2] = 2.0*numpy.sum(y1[indices1]/lam_1[indices1]**2) \
								  + 2.0*numpy.sum(y2[indices2]/lam_2[indices2]**2)

			hess_res[3*iq+1,3*iq] = hess_res[3*iq,3*iq+1]
			hess_res[3*iq+1,3*iq+1] = 2.0*numpy.sum((x1[indices1]/lam_1[indices1])**2) \
									+ 2.0*numpy.sum((x2[indices2]/lam_2[indices2])**2)
			hess_res[3*iq+1,3*iq+2] = 2.0*numpy.sum(x1[indices1]*y1[indices1]/lam_1[indices1]**2) \
									+ 2.0*numpy.sum(x2[indices2]*y2[indices2]/lam_2[indices2]**2)

			hess_res[3*iq+2,3*iq] = hess_res[3*iq,3*iq+2]
			hess_res[3*iq+2,3*iq+1] = hess_res[3*iq+1,3*iq+2]
			hess_res[3*iq+2,3*iq+2] = 2.0*numpy.sum(y1[indices1]**2/lam_1[indices1]**2) \
									+ 2.0*numpy.sum(y2[indices2]**2/lam_2[indices2]**2)

			hess_res[12+3*iq,12+3*iq] = hess_res[3*iq,3*iq]
			hess_res[12+3*iq,12+3*iq+1] = hess_res[3*iq,3*iq+1]
			hess_res[12+3*iq,12+3*iq+2] = hess_res[3*iq,3*iq+2]
			hess_res[12+3*iq+1,12+3*iq] = hess_res[3*iq+1,3*iq]
			hess_res[12+3*iq+1,12+3*iq+1] = hess_res[3*iq+1,3*iq+1]
			hess_res[12+3*iq+1,12+3*iq+2] = hess_res[3*iq+1,3*iq+2]
			hess_res[12+3*iq+2,12+3*iq] = hess_res[3*iq+2,3*iq]
			hess_res[12+3*iq+2,12+3*iq+1] = hess_res[3*iq+2,3*iq+1]
			hess_res[12+3*iq+2,12+3*iq+2] = hess_res[3*iq+2,3*iq+2]

		for iq in range(4):

			indices1 = numpy.argwhere(q1==iq)
			indices2 = numpy.argwhere(q2==iq)

			hess_res[3*iq,24+2*iq] = - 4*numpy.sum(x1[indices1]*new_pts[indices1,0]/lam_1[indices1]**2) \
							       	 - 4*numpy.sum(x2[indices2]*new_pts[indices2,3]/lam_2[indices2]**2) \
									 + 2*numpy.sum(x1[indices1]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	 + 2*numpy.sum(x2[indices2]*new_pts[indices2,0]/lam_2[indices2]**2)

			hess_res[3*iq,24+2*iq+1] = - 4*numpy.sum(y1[indices1]*new_pts[indices1,0]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(y2[indices2]*new_pts[indices2,3]/lam_2[indices2]**2) \
									   + 2*numpy.sum(y1[indices1]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(y2[indices2]*new_pts[indices2,0]/lam_2[indices2]**2)

			hess_res[3*iq+1,24+2*iq] = - 4*numpy.sum(x1[indices1]**2*new_pts[indices1,0]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(x2[indices2]**2*new_pts[indices2,3]/lam_2[indices2]**2) \
									   + 2*numpy.sum(x1[indices1]**2*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(x2[indices2]**2*new_pts[indices2,0]/lam_2[indices2]**2)

			hess_res[3*iq+1,24+2*iq+1] = - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,0]/lam_1[indices1]**2) \
							       	     - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,3]/lam_2[indices2]**2) \
									     + 2*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	     + 2*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,0]/lam_2[indices2]**2)

			hess_res[3*iq+2,24+2*iq] = - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,0]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,3]/lam_2[indices2]**2) \
									   + 2*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,0]/lam_2[indices2]**2)

			hess_res[3*iq+2,24+2*iq+1] = - 4*numpy.sum(y1[indices1]**2*new_pts[indices1,0]/lam_1[indices1]**2) \
							        	 - 4*numpy.sum(y2[indices2]**2*new_pts[indices2,3]/lam_2[indices2]**2) \
									     + 2*numpy.sum(y1[indices1]**2*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	     + 2*numpy.sum(y2[indices2]**2*new_pts[indices2,0]/lam_2[indices2]**2)


		for iq in range(4):

			indices1 = numpy.argwhere(q1==iq)
			indices2 = numpy.argwhere(q2==iq)

			hess_res[12+3*iq,24+2*iq] = - 4*numpy.sum(x1[indices1]*new_pts[indices1,1]/lam_1[indices1]**2) \
							       	 - 4*numpy.sum(x2[indices2]*new_pts[indices2,4]/lam_2[indices2]**2) \
									 + 2*numpy.sum(x1[indices1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	 + 2*numpy.sum(x2[indices2]*new_pts[indices2,1]/lam_2[indices2]**2)

			hess_res[12+3*iq,24+2*iq+1] = - 4*numpy.sum(y1[indices1]*new_pts[indices1,1]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(y2[indices2]*new_pts[indices2,4]/lam_2[indices2]**2) \
									   + 2*numpy.sum(y1[indices1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(y2[indices2]*new_pts[indices2,1]/lam_2[indices2]**2)

			hess_res[12+3*iq+1,24+2*iq] = - 4*numpy.sum(x1[indices1]**2*new_pts[indices1,1]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(x2[indices2]**2*new_pts[indices2,4]/lam_2[indices2]**2) \
									   + 2*numpy.sum(x1[indices1]**2*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(x2[indices2]**2*new_pts[indices2,1]/lam_2[indices2]**2)

			hess_res[12+3*iq+1,24+2*iq+1] = - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,1]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,4]/lam_2[indices2]**2) \
									   + 2*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,1]/lam_2[indices2]**2)

			hess_res[12+3*iq+2,24+2*iq] = - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,1]/lam_1[indices1]**2) \
							       	   - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,4]/lam_2[indices2]**2) \
									   + 2*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	   + 2*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,1]/lam_2[indices2]**2)

			hess_res[12+3*iq+2,24+2*iq+1] = - 4*numpy.sum(y1[indices1]**2*new_pts[indices1,1]/lam_1[indices1]**2) \
							        	 - 4*numpy.sum(y2[indices2]**2*new_pts[indices2,4]/lam_2[indices2]**2) \
									     + 2*numpy.sum(y1[indices1]**2*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	     + 2*numpy.sum(y2[indices2]**2*new_pts[indices2,1]/lam_2[indices2]**2)


		for iq in range(4):

			indices1 = numpy.argwhere(q1==iq)
			indices2 = numpy.argwhere(q2==iq)

			hess_res[24+2*iq,24+2*iq] = + 6*numpy.sum(x1[indices1]**2*new_pts[indices1,0]**2/lam_1[indices1]**2) \
							     	    + 6*numpy.sum(x2[indices2]**2*new_pts[indices2,3]**2/lam_2[indices2]**2) \
									    - 4*numpy.sum(x1[indices1]**2*new_pts[indices1,0]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	    - 4*numpy.sum(x2[indices2]**2*new_pts[indices2,0]*new_pts[indices2,3]/lam_2[indices2]**2) \
							       	    + 6*numpy.sum(x1[indices1]**2*new_pts[indices1,1]**2/lam_1[indices1]**2) \
							     	    + 6*numpy.sum(x2[indices2]**2*new_pts[indices2,4]**2/lam_2[indices2]**2) \
									    - 4*numpy.sum(x1[indices1]**2*new_pts[indices1,1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	    - 4*numpy.sum(x2[indices2]**2*new_pts[indices2,1]*new_pts[indices2,4]/lam_2[indices2]**2)

			hess_res[24+2*iq,24+2*iq+1] = + 6*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,0]**2/lam_1[indices1]**2) \
							     	      + 6*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,3]**2/lam_2[indices2]**2) \
									      - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,0]*new_pts[indices1,3]/lam_1[indices1]**2) \
							       	      - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,0]*new_pts[indices2,3]/lam_2[indices2]**2) \
							       	      + 6*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,1]**2/lam_1[indices1]**2) \
							     	      + 6*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,4]**2/lam_2[indices2]**2) \
									      - 4*numpy.sum(x1[indices1]*y1[indices1]*new_pts[indices1,1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							        	  - 4*numpy.sum(x2[indices2]*y2[indices2]*new_pts[indices2,1]*new_pts[indices2,4]/lam_2[indices2]**2)

			hess_res[24+2*iq+1,24+2*iq] = hess_res[24+2*iq,24+2*iq+1]

			hess_res[24+2*iq+1,24+2*iq+1] = + 6*numpy.sum(y1[indices1]**2*new_pts[indices1,0]**2/lam_1[indices1]**2) \
							     	        + 6*numpy.sum(y2[indices2]**2*new_pts[indices2,3]**2/lam_2[indices2]**2) \
									        - 4*numpy.sum(y1[indices1]**2*new_pts[indices1,0]*new_pts[indices1,3]/lam_1[indices1]**2) \
							         	    - 4*numpy.sum(y2[indices2]**2*new_pts[indices2,0]*new_pts[indices2,3]/lam_2[indices2]**2) \
							        	    + 6*numpy.sum(y1[indices1]**2*new_pts[indices1,1]**2/lam_1[indices1]**2) \
							     	        + 6*numpy.sum(y2[indices2]**2*new_pts[indices2,4]**2/lam_2[indices2]**2) \
									        - 4*numpy.sum(y1[indices1]**2*new_pts[indices1,1]*new_pts[indices1,4]/lam_1[indices1]**2) \
							       	        - 4*numpy.sum(y2[indices2]**2*new_pts[indices2,1]*new_pts[indices2,4]/lam_2[indices2]**2)

		for iq1 in range(4):

			for iq2 in range(4):

				if iq1==iq2: continue

				indices12 = numpy.argwhere((q1==iq1) & (q2==iq2))
				indices21 = numpy.argwhere((q1==iq2) & (q2==iq1))

				hess_res[3*iq1,3*iq2] = - 2.0*numpy.sum(1.0/lam_1[indices12]/lam_2[indices12]) \
										- 2.0*numpy.sum(1.0/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1,3*iq2+1] = - 2.0*numpy.sum(x2[indices12]/lam_1[indices12]/lam_2[indices12]) \
										  - 2.0*numpy.sum(x1[indices21]/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1,3*iq2+2] = - 2.0*numpy.sum(y2[indices12]/lam_1[indices12]/lam_2[indices12]) \
										  - 2.0*numpy.sum(y1[indices21]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1+1,3*iq2] = - 2.0*numpy.sum(x1[indices12]/lam_1[indices12]/lam_2[indices12]) \
										  - 2.0*numpy.sum(x2[indices21]/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1+1,3*iq2+1] = - 2.0*numpy.sum(x1[indices12]*x2[indices12]/lam_1[indices12]/lam_2[indices12]) \
											- 2.0*numpy.sum(x2[indices21]*x1[indices21]/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1+1,3*iq2+2] = - 2.0*numpy.sum(x1[indices12]*y2[indices12]/lam_1[indices12]/lam_2[indices12]) \
				 							- 2.0*numpy.sum(x2[indices21]*y1[indices21]/lam_1[indices21]/lam_2[indices21])


				hess_res[3*iq1+2,3*iq2] = - 2.0*numpy.sum(y1[indices12]/lam_1[indices12]/lam_2[indices12]) \
				 						  - 2.0*numpy.sum(y2[indices21]/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1+2,3*iq2+1] = - 2.0*numpy.sum(y1[indices12]*x2[indices12]/lam_1[indices12]/lam_2[indices12]) \
											- 2.0*numpy.sum(y2[indices21]*x1[indices21]/lam_1[indices21]/lam_2[indices21])
				hess_res[3*iq1+2,3*iq2+2] = - 2.0*numpy.sum(y1[indices12]*y2[indices12]/lam_1[indices12]/lam_2[indices12]) \
											- 2.0*numpy.sum(y2[indices21]*y1[indices21]/lam_1[indices21]/lam_2[indices21])


		for iq1 in range(4):

			for iq2 in range(4):

				if iq1==iq2: continue

				indices12 = numpy.argwhere((q1==iq1) & (q2==iq2))
				indices21 = numpy.argwhere((q1==iq2) & (q2==iq1))

				hess_res[3*iq1,24+2*iq2] =  2*numpy.sum(x2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								          + 2*numpy.sum(x1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1,24+2*iq2+1] =  2*numpy.sum(y2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								            + 2*numpy.sum(y1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1+1,24+2*iq2] =  2*numpy.sum(x1[indices12]*x2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								            + 2*numpy.sum(x2[indices21]*x1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1+1,24+2*iq2+1] =  2*numpy.sum(x1[indices12]*y2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								              + 2*numpy.sum(x2[indices21]*y1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1+2,24+2*iq2] =  2*numpy.sum(y1[indices12]*x2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								            + 2*numpy.sum(y2[indices21]*x1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])

				hess_res[3*iq1+2,24+2*iq2+1] =  2*numpy.sum(y1[indices12]*y2[indices12]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								              + 2*numpy.sum(y2[indices21]*y1[indices21]*new_pts[indices21,0]/lam_1[indices21]/lam_2[indices21])


		for iq1 in range(4):

			for iq2 in range(4):

				if iq1==iq2: continue

				indices12 = numpy.argwhere((q1==iq1) & (q2==iq2))
				indices21 = numpy.argwhere((q1==iq2) & (q2==iq1))

				hess_res[12+3*iq1,24+2*iq2] =  2*numpy.sum(x1[indices12]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								          + 2*numpy.sum(x2[indices21]*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])

				hess_res[12+3*iq1,24+2*iq2+1] =  2*numpy.sum(y1[indices12]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								            + 2*numpy.sum(y2[indices21]*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])


				hess_res[12+3*iq1+1,24+2*iq2] =  2*numpy.sum(x1[indices12]**2*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								          + 2*numpy.sum(x2[indices21]**2*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])

				hess_res[12+3*iq1+1,24+2*iq2+1] =  2*numpy.sum(x1[indices12]*y1[indices12]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								              + 2*numpy.sum(x2[indices21]*y2[indices21]*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])

				hess_res[12+3*iq1+2,24+2*iq2] =  2*numpy.sum(x1[indices12]*y1[indices12]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								            + 2*numpy.sum(x2[indices21]*y2[indices21]*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])

				hess_res[12+3*iq1+2,24+2*iq2+1] =  2*numpy.sum(y1[indices12]**2*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								              + 2*numpy.sum(y2[indices21]**2*new_pts[indices21,1]/lam_1[indices21]/lam_2[indices21])

		for iq1 in range(4):

			for iq2 in range(4):

				if iq1==iq2: continue

				indices12 = numpy.argwhere((q1==iq1) & (q2==iq2))
				indices21 = numpy.argwhere((q1==iq2) & (q2==iq1))

				hess_res[24+2*iq1,24+2*iq2] = - 2*numpy.sum(x1[indices12]*x2[indices12]*new_pts[indices12,0]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								              - 2*numpy.sum(x1[indices21]*x2[indices21]*new_pts[indices21,0]*new_pts[indices21,3]/lam_1[indices21]/lam_2[indices21]) \
								              - 2*numpy.sum(x1[indices12]*x2[indices12]*new_pts[indices12,1]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								              - 2*numpy.sum(x1[indices21]*x2[indices21]*new_pts[indices21,1]*new_pts[indices21,4]/lam_1[indices21]/lam_2[indices21])

				hess_res[24+2*iq1,24+2*iq2+1] = - 2*numpy.sum(x1[indices12]*y2[indices12]*new_pts[indices12,0]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								                - 2*numpy.sum(x1[indices21]*y2[indices21]*new_pts[indices21,0]*new_pts[indices21,3]/lam_1[indices21]/lam_2[indices21]) \
								                - 2*numpy.sum(x1[indices12]*y2[indices12]*new_pts[indices12,1]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								                - 2*numpy.sum(x1[indices21]*y2[indices21]*new_pts[indices21,1]*new_pts[indices21,4]/lam_1[indices21]/lam_2[indices21])

				hess_res[24+2*iq1+1,24+2*iq2] = hess_res[24+2*iq1,24+2*iq2+1]

				hess_res[24+2*iq1+1,24+2*iq2+1] = - 2*numpy.sum(y1[indices12]*y2[indices12]*new_pts[indices12,0]*new_pts[indices12,3]/lam_1[indices12]/lam_2[indices12]) \
								                  - 2*numpy.sum(y1[indices21]*y2[indices21]*new_pts[indices21,0]*new_pts[indices21,3]/lam_1[indices21]/lam_2[indices21]) \
								             	  - 2*numpy.sum(y1[indices12]*y2[indices12]*new_pts[indices12,1]*new_pts[indices12,4]/lam_1[indices12]/lam_2[indices12]) \
								                  - 2*numpy.sum(y1[indices21]*y2[indices21]*new_pts[indices21,1]*new_pts[indices21,4]/lam_1[indices21]/lam_2[indices21])


		hess_res[12:24,12:24] = hess_res[:12,:12]  #False

		hess_res[24:,:24] = hess_res[:24,24:].T

		for iq in q_fix:

			if iq%2==0:

				X = Tx[3*iq] + Tx[3*iq+1]*ext[iq,0] + Tx[3*iq+2]*ext[iq,1]
				lam = 1.0 + lambda_[2*iq]*ext[iq,0] + lambda_[2*iq+1]*ext[iq,1]
				lx = X/lam - ext[iq,0]

				hess_res[3*iq,3*iq] += 2.0/lam**2
				hess_res[3*iq,3*iq+1] += 2.0*ext[iq,0]/lam**2
				hess_res[3*iq,3*iq+2] += 2.0*ext[iq,1]/lam**2
				hess_res[3*iq+1,3*iq] += 2.0*ext[iq,0]/lam**2
				hess_res[3*iq+1,3*iq+1] += 2.0*ext[iq,0]**2/lam**2
				hess_res[3*iq+1,3*iq+2] += 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[3*iq+2,3*iq] += 2.0*ext[iq,1]/lam**2
				hess_res[3*iq+2,3*iq+1] += 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[3*iq+2,3*iq+2] += 2.0*ext[iq,1]**2/lam**2

				hess_res[3*iq,24+2*iq] += - 4.0*X*ext[iq,0]/lam**3 + 2.0*ext[iq,0]**2/lam**2
				hess_res[3*iq,24+2*iq+1] += - 4.0*X*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[3*iq+1,24+2*iq] += - 4.0*X*ext[iq,0]**2/lam**3 + 2.0*ext[iq,0]**3/lam**2
				hess_res[3*iq+1,24+2*iq+1] += - 4.0*X*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[3*iq+2,24+2*iq] += - 4.0*X*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[3*iq+2,24+2*iq+1] += - 4.0*X*ext[iq,1]**2/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2

				hess_res[24+2*iq,3*iq] += - 4.0*X*ext[iq,0]/lam**3 + 2.0*ext[iq,0]**2/lam**2
				hess_res[24+2*iq+1,3*iq] += - 4.0*X*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[24+2*iq,3*iq+1] += - 4.0*X*ext[iq,0]**2/lam**3 + 2.0*ext[iq,0]**3/lam**2
				hess_res[24+2*iq+1,3*iq+1] += - 4.0*X*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[24+2*iq,3*iq+2] += - 4.0*X*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[24+2*iq+1,3*iq+2] += - 4.0*X*ext[iq,1]**2/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2

				hess_res[24+2*iq,24+2*iq] += 6.0*X**2*ext[iq,0]**2/lam**4 - 4.0*X*ext[iq,0]**3/lam**3
				hess_res[24+2*iq,24+2*iq+1] += 6.0*X**2*ext[iq,0]*ext[iq,1]/lam**4 - 4.0*X*ext[iq,0]**2*ext[iq,1]/lam**3
				hess_res[24+2*iq+1,24+2*iq] += 6.0*X**2*ext[iq,0]*ext[iq,1]/lam**4 - 4.0*X*ext[iq,0]**2*ext[iq,1]/lam**3
				hess_res[24+2*iq+1,24+2*iq+1] += 6.0*X**2*ext[iq,1]**2/lam**4 - 4.0*X*ext[iq,0]*ext[iq,1]**2/lam**3


			if iq%2==1:

				Y = Ty[3*iq] + Ty[3*iq+1]*ext[iq,0] + Ty[3*iq+2]*ext[iq,1]
				lam = 1.0 + lambda_[2*iq]*ext[iq,0] + lambda_[2*iq+1]*ext[iq,1]
				ly = Y/lam - ext[iq,1]

				hess_res[12+3*iq,12+3*iq] += 2.0/lam**2
				hess_res[12+3*iq,12+3*iq+1] += 2.0*ext[iq,0]/lam**2
				hess_res[12+3*iq,12+3*iq+2] += 2.0*ext[iq,1]/lam**2
				hess_res[12+3*iq+1,12+3*iq] += 2.0*ext[iq,0]/lam**2
				hess_res[12+3*iq+1,12+3*iq+1] += 2.0*ext[iq,0]**2/lam**2
				hess_res[12+3*iq+1,12+3*iq+2] += 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[12+3*iq+2,12+3*iq] += 2.0*ext[iq,1]/lam**2
				hess_res[12+3*iq+2,12+3*iq+1] += 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[12+3*iq+2,12+3*iq+2] += 2.0*ext[iq,1]**2/lam**2

				hess_res[12+3*iq,24+2*iq] += - 4.0*Y*ext[iq,0]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[12+3*iq,24+2*iq+1] += - 4.0*Y*ext[iq,1]/lam**3 + 2.0*ext[iq,1]**2/lam**2
				hess_res[12+3*iq+1,24+2*iq] += - 4.0*Y*ext[iq,0]**2/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[12+3*iq+1,24+2*iq+1] += - 4.0*Y*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2
				hess_res[12+3*iq+2,24+2*iq] += - 4.0*Y*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2
				hess_res[12+3*iq+2,24+2*iq+1] += - 4.0*Y*ext[iq,1]**2/lam**3 + 2.0*ext[iq,1]**3/lam**2

				hess_res[24+2*iq,12+3*iq] += - 4.0*Y*ext[iq,0]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]/lam**2
				hess_res[24+2*iq+1,12+3*iq] += - 4.0*Y*ext[iq,1]/lam**3 + 2.0*ext[iq,1]**2/lam**2
				hess_res[24+2*iq,12+3*iq+1] += - 4.0*Y*ext[iq,0]**2/lam**3 + 2.0*ext[iq,0]**2*ext[iq,1]/lam**2
				hess_res[24+2*iq+1,12+3*iq+1] += - 4.0*Y*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2
				hess_res[24+2*iq,12+3*iq+2] += - 4.0*Y*ext[iq,0]*ext[iq,1]/lam**3 + 2.0*ext[iq,0]*ext[iq,1]**2/lam**2
				hess_res[24+2*iq+1,12+3*iq+2] += - 4.0*Y*ext[iq,1]**2/lam**3 + 2.0*ext[iq,1]**3/lam**2

				hess_res[24+2*iq,24+2*iq] += 6.0*Y**2*ext[iq,0]**2/lam**4 - 4.0*Y*ext[iq,0]**2*ext[iq,1]/lam**3
				hess_res[24+2*iq,24+2*iq+1] += 6.0*Y**2*ext[iq,0]*ext[iq,1]/lam**4 - 4.0*Y*ext[iq,0]*ext[iq,1]**2/lam**3
				hess_res[24+2*iq+1,24+2*iq] += 6.0*Y**2*ext[iq,0]*ext[iq,1]/lam**4 - 4.0*Y*ext[iq,0]*ext[iq,1]**2/lam**3
				hess_res[24+2*iq+1,24+2*iq+1] += 6.0*Y**2*ext[iq,1]**2/lam**4 - 4.0*Y*ext[iq,1]**3/lam**3

		return hess_res

	def res(T,q_fix,ext):

		r = residual(T,q_fix,ext)

		return r[1][2]


	# Evaluate the possible anchor points

	(Tx_,Ty_,new_pts_,res_) = getPureTranslations(pts)

	xy = numpy.squeeze((new_pts_[:,:2] + new_pts_[:,3:5])/2.0)

	ext1 = xy[numpy.argmin(xy[:,0])]
	ext2 = xy[numpy.argmax(xy[:,1])]
	ext3 = xy[numpy.argmax(xy[:,0])]
	ext4 = xy[numpy.argmin(xy[:,1])]

	ext = numpy.row_stack((ext1,ext2,ext3,ext4))

	xmin0 = numpy.max(Tx_[0], Tx_[9])
	ymin0 = numpy.max(Ty_[6], Ty_[9])

	xlength = numpy.min(2*nx + WIDTH - Tx_[0] + Tx_[3], 2*nx - WIDTH + Tx_[6] - Tx_[9])
	ylength = numpy.min(2*ny + WIDTH - Ty_[9] + Ty_[0], 2*ny - WIDTH + Ty_[3] - Ty_[6])

	frame = numpy.array([[ xmin0, xmin0 + xlength - 1 ],[ ymin0, ymin0 + ylength - 1 ]])

	# Initial estimate

	T0 = numpy.zeros((32))

	T0[0:12] = Tx_[:]
	T0[12:24] = Ty_[:]

#	X = numpy.random.random((32))
#
#	util.tester.test_derivatives(res,der_residual,X,args=(qfix,ext),delta_h=1e-5,d2fdX2=hess_res)
#
#	raise SystemExit

	if use_ncg:
		T0 = scipy.optimize.fmin_ncg(res, T0, der_residual, args = (qfix, ext), fhess=hess_res, avextol=1.e-16, maxiter=ITERATION_NUMBER )
	else:
		T0 = scipy.optimize.fmin(res, T0, args = (qfix, ext), xtol=1e-6, ftol=1e-6, maxiter=ITERATION_NUMBER, maxfun=ITERATION_NUMBER)

	Tx,Ty = T0[:12],T0[12:24]
	lambda_x,lambda_y = T0[24:32],T0[24:32]

	# Reshuffle Tx and Ty

	Tx = numpy.concatenate((Tx,lambda_x))
	Ty = numpy.concatenate((Ty,lambda_y))

	(new_pts, res, resx, resy) = residual(T0, qfix, ext)

	log.info('Residual:')
	for index,pt in enumerate(zip(resx,resy,new_pts[:,0],new_pts[:,1])):
		log.info('Stitch #%i:	Res. in x: %- 10.3f	Res in y: %- 10.3f	(x,y)=(%5.1f,%5.1f)' %(index,pt[0],pt[1],pt[2],pt[3]))

	(Tx,Ty,frame) = rescale(Tx,Ty,frame,nx,ny)

	return (Tx,Ty,new_pts,res,frame)


def getProjectiveTransformations_2(pts, nx, ny, qfix=[0,1]):

	global iter
	iter = 0

	def residual(T,q_fix,ext):

		Tx = T[0:12]
		Ty = T[12:24]
		lambda_x = T[24:32]
		lambda_y = T[32:40]

		x1,y1,q1 = pts[:,0],pts[:,1],pts[:,2].astype('int')
		x2,y2,q2 = pts[:,3],pts[:,4],pts[:,5].astype('int')

		new_pts = numpy.zeros_like(pts)

		lam_1x = 1.0 + lambda_x[2*q1]*x1 + lambda_x[2*q1+1]*y1
		lam_1y = 1.0 + lambda_y[2*q1]*x1 + lambda_y[2*q1+1]*y1

		lam_2x = 1.0 + lambda_x[2*q2]*x2 + lambda_x[2*q2+1]*y2
		lam_2y = 1.0 + lambda_y[2*q2]*x2 + lambda_y[2*q2+1]*y2

		new_pts[:,0] = (Tx[3*q1] + Tx[3*q1+1]*x1 + Tx[3*q1+2]*y1)/lam_1x
		new_pts[:,1] = (Ty[3*q1] + Ty[3*q1+1]*x1 + Ty[3*q1+2]*y1)/lam_1y

		new_pts[:,3] = (Tx[3*q2] + Tx[3*q2+1]*x2 + Tx[3*q2+2]*y2)/lam_2x
		new_pts[:,4] = (Ty[3*q2] + Ty[3*q2+1]*x2 + Ty[3*q2+2]*y2)/lam_2y

		resx = new_pts[:,0] - new_pts[:,3]
		resy = new_pts[:,1] - new_pts[:,4]

		res2 = resx**2 + resy**2
		res = numpy.sqrt(res2)

		global iter

		line_output = ' #%-8i	Mean Res.: %-8.3f	Max. Res.: %-8.3f	Full Res. Sqr.: %-8.3f' %(iter,numpy.mean(res), numpy.max(res), numpy.sum(res2))

		if iter%NPRINT==0:
			print "\r%s" %line_output
		iter += 1

		res2 = numpy.sum(res2)

		for iq in q_fix:

			lx = (Tx[3*iq] + Tx[3*iq+1]*ext[iq,0] + Tx[3*iq+2]*ext[iq,1])/(1.0 + lambda_x[2*iq]*ext[iq,0] + lambda_x[2*iq+1]*ext[iq,1]) - ext[iq,0]
			ly = (Ty[3*iq] + Ty[3*iq+1]*ext[iq,0] + Ty[3*iq+2]*ext[iq,1])/(1.0 + lambda_y[2*iq]*ext[iq,0] + lambda_y[2*iq+1]*ext[iq,1]) - ext[iq,1]

			if iq%2==0: res2 += lx**2
			if iq%2==1: res2 += ly**2

		res = numpy.array([numpy.min(res), numpy.max(res), numpy.sum(res2), numpy.mean(res)])

		return (new_pts,res,resx,resy)


	def res(T,q_fix,ext):

		r = residual(T,q_fix,ext)

		return r[1][2]

	# Evaluate the possible anchor points

	(Tx_,Ty_,new_pts_,res_) = getPureTranslations(pts)

	xy = numpy.squeeze((new_pts_[:,:2] + new_pts_[:,3:5])/2.0)

	ext1 = xy[numpy.argmin(xy[:,0])]
	ext2 = xy[numpy.argmax(xy[:,1])]
	ext3 = xy[numpy.argmax(xy[:,0])]
	ext4 = xy[numpy.argmin(xy[:,1])]

	ext = numpy.row_stack((ext1,ext2,ext3,ext4))

	xmin0 = numpy.max(Tx_[0], Tx_[9])
	ymin0 = numpy.max(Ty_[6], Ty_[9])

	xlength = numpy.min(2*nx + WIDTH - Tx_[0] + Tx_[3], 2*nx - WIDTH + Tx_[6] - Tx_[9])
	ylength = numpy.min(2*ny + WIDTH - Ty_[9] + Ty_[0], 2*ny - WIDTH + Ty_[3] - Ty_[6])

	frame = numpy.array([[ xmin0, xmin0 + xlength - 1 ],[ ymin0, ymin0 + ylength - 1 ]])

	# Initial estimate

	T0 = numpy.zeros((40))	# First estimate

	T0[:12] = Tx_[:]
	T0[12:24] = Ty_[:]

	T0 = scipy.optimize.fmin(res, T0, args=(qfix,ext,), xtol=1e-6, ftol=1e-6, maxiter=ITERATION_NUMBER, maxfun=ITERATION_NUMBER)

	Tx,Ty = T0[:12],T0[12:24]
	lambda_x,lambda_y = T0[24:32],T0[32:40]

	(new_pts,res, resx, resy) = residual(T0,qfix,ext)

	log.info('Residual:')
	for index,pt in enumerate(zip(resx,resy,new_pts[:,0],new_pts[:,1])):
		log.info('Stitch #%i:	Res. in x: %- 10.3f	Res in y: %- 10.3f	(x,y)=(%5.1f,%5.1f)' %(index,pt[0],pt[1],pt[2],pt[3]))

	# Reshuffle Tx and Ty

	Tx = numpy.concatenate((Tx,lambda_x))
	Ty = numpy.concatenate((Ty,lambda_y))

	(Tx,Ty,frame) = rescale(Tx,Ty,frame,nx,ny)

	return (Tx,Ty,new_pts,res,frame)


def rescale( Tx,Ty,frame,nx,ny):

	Tx_,Ty_ = Tx[:12],Ty[:12]
	lambda_x_,lambda_y_ = Tx[12:],Ty[12:]

	log.info('Rescaling the transformation!')

	xmin0 = frame[0,0]
	ymin0 = frame[1,0]
	xlength = frame[0,1] - frame[0,0]
	ylength = frame[1,1] - frame[1,0]

	base_frame = numpy.array([[1,ny],[nx,ny],[nx,1],[1,1]])

	frame_old = numpy.row_stack(( \
								  base_frame + [0,ny+WIDTH], \
								  base_frame + [nx+WIDTH,ny+WIDTH], \
								  base_frame + [nx+WIDTH,0], \
								  base_frame + [0,0] \
								  ))

	xmin_old = (frame_old[0,0] + frame_old[15,0])/2.0
	xmax_old = (frame_old[5,0] + frame_old[10,0])/2.0
	ymin_old = (frame_old[10,1] + frame_old[14,1])/2.0
	ymax_old = (frame_old[0,1] + frame_old[5,1])/2.0

	# Evaluate the image of the frame (for the 4 quadrants)

	frame_new = numpy.empty_like(frame_old)

	for iq in range(4):
		x = frame_old[4*iq:4*iq+4,0]
		y = frame_old[4*iq:4*iq+4,1]
		lamx = 1.0 + lambda_x_[2*iq]*x + lambda_x_[2*iq+1]*y
		lamy = 1.0 + lambda_y_[2*iq]*x + lambda_y_[2*iq+1]*y
		frame_new[4*iq:4*iq+4,0] = (Tx_[3*iq] + Tx_[3*iq+1]*x + Tx_[3*iq+2]*y)/lamx
		frame_new[4*iq:4*iq+4,1] = (Ty_[3*iq] + Ty_[3*iq+1]*x + Ty_[3*iq+2]*y)/lamy

	xmin_new = (frame_new[0,0] + frame_new[15,0])/2.0
	xmax_new = (frame_new[5,0] + frame_new[10,0])/2.0
	ymin_new = (frame_new[10,1] + frame_new[14,1])/2.0
	ymax_new = (frame_new[0,1] + frame_new[5,1])/2.0

	log.info('Frame for the original transfromation:')
	for iq in range(4):
		log.info('Quadrant %i:' %iq)
		log.info('      %s' %(frame_old[4*iq:4*iq+4].tolist()))
		log.info('   -> %s' %(frame_new[4*iq:4*iq+4].tolist()))

	log.info('(xmin,xmax): (%.1f,%.1f)->(%.1f,%.1f)' %(xmin_old,xmax_old,xmin_new,xmax_new))
	log.info('(ymin,ymax): (%.1f,%.1f)->(%.1f,%.1f)' %(ymin_old,ymax_old,ymin_new,ymax_new))

	# Rescaling Parameters

	xmin_old = xmin0
	ymin_old = ymin0

	lam_x = xlength/(xmax_new-xmin_new)
	lam_y = ylength/(ymax_new-ymin_new)

	log.info('Scaling Parameters for the new transformation in (x,y) = (%f,%f)' %(lam_x,lam_y))

	# Rescale the transformation

	for iq in range(4):

		Tx_[3*iq] = lam_x*(Tx_[3*iq] - xmin_new) + xmin_old
		Tx_[3*iq+1] = lam_x*(Tx_[3*iq+1] - xmin_new*lambda_x_[2*iq]) + xmin_old*lambda_x_[2*iq]
		Tx_[3*iq+2] = lam_x*(Tx_[3*iq+2] - xmin_new*lambda_x_[2*iq+1]) + xmin_old*lambda_x_[2*iq+1]

		Ty_[3*iq] = lam_y*(Ty_[3*iq] - ymin_new) + ymin_old
		Ty_[3*iq+1] = lam_y*(Ty_[3*iq+1] - ymin_new*lambda_y_[2*iq]) + ymin_old*lambda_y_[2*iq]
		Ty_[3*iq+2] = lam_y*(Ty_[3*iq+2] - ymin_new*lambda_y_[2*iq+1]) + ymin_old*lambda_y_[2*iq+1]

	# Evaluate the new frame (for the 4 quadrants). This is the image of the initial frame
	# with the rescaled transfromation

	for iq in range(4):
		x = frame_old[4*iq:4*iq+4,0]
		y = frame_old[4*iq:4*iq+4,1]
		lamx = 1.0 + lambda_x_[2*iq]*x + lambda_x_[2*iq+1]*y
		lamy = 1.0 + lambda_y_[2*iq]*x + lambda_y_[2*iq+1]*y
		frame_new[4*iq:4*iq+4,0] = (Tx_[3*iq] + Tx_[3*iq+1]*x + Tx_[3*iq+2]*y)/lamx
		frame_new[4*iq:4*iq+4,1] = (Ty_[3*iq] + Ty_[3*iq+1]*x + Ty_[3*iq+2]*y)/lamy

	xmin = numpy.max([frame_new[0,0],frame_new[3,0],frame_new[12,0],frame_new[15,0]])
	xmax = numpy.min([frame_new[5,0],frame_new[6,0],frame_new[9,0],frame_new[10,0]])

	ymin = numpy.max([frame_new[10,1],frame_new[11,1],frame_new[14,1],frame_new[15,1]])
	ymax = numpy.min([frame_new[0,1],frame_new[1,1],frame_new[4,1],frame_new[5,1]])

	xmin = numpy.ceil(xmin+1)
	xmax = numpy.floor(xmax-1)

	ymin = numpy.ceil(ymin+1)
	ymax = numpy.floor(ymax-1)

	frame = numpy.array([[xmin,xmax],[ymin,ymax]])

	log.info('Frame for the rescaled transfromation:')
	for iq in range(4):
		log.info('Quadrant %i:' %iq)
		log.info('      %s' %(frame_old[4*iq:4*iq+4].tolist()))
		log.info('   -> %s' %(frame_new[4*iq:4*iq+4].tolist()))

	log.info('(xmin,xmax): (%.1f,%.1f)->(%.1f,%.1f)' %(xmin_old,xmax_old,xmin_new,xmax_new))
	log.info('(ymin,ymax): (%.1f,%.1f)->(%.1f,%.1f)' %(ymin_old,ymax_old,ymin_new,ymax_new))

	return (Tx, Ty, frame)


def inverseTransformation(Tx,Ty,nx,ny,mode):
	'''Calculate the transformtions. Inverse the map.'''

	if mode==TRANSLATIONS or mode==AFFINE_TRANSFORMATIONS:

		R = numpy.zeros((2,2))

		for i in range(4):

			R[0,:] = Tx[3*i+1:3*i+3]
			R[1,:] = Ty[3*i+1:3*i+3]

			R = numpy.linalg.pinv(R)

			Tx[3*i] = - R[0,0]*Tx[3*i] - R[0,1]*Ty[3*i]
			Tx[3*i+1] = R[0,0]
			Tx[3*i+2] = R[0,1]

			Ty[3*i] = - R[1,0]*Tx[3*i] - R[1,1]*Ty[3*i]
			Ty[3*i+1] = R[1,0]
			Ty[3*i+2] = R[1,1]

		Tx[3] = Tx[3] - nx - WIDTH
		Tx[6] = Tx[6] - nx - WIDTH

		Ty[0] = Ty[0] - ny - WIDTH
		Ty[3] = Ty[3] - ny - WIDTH

		transformations = numpy.row_stack((Tx,Ty))

		return transformations

	elif mode==PROJECTIVE_TRANSFORMATIONS_1 or mode==PROJECTIVE_TRANSFORMATIONS_2:

		Tx_ = Tx[0:12]
		lambda_x_ = Tx[12:]

		Ty_ = Ty[0:12]
		lambda_y_ = Ty[12:]

		Tx_inv = numpy.zeros((16))
		Ty_inv = numpy.zeros((16))
		lam_inv = numpy.zeros((16))	# should be 12 normally

		for i in range(4):

			tx,ty = Tx_[3*i],Ty_[3*i]

			lam_00,lam_01,lam_10,lam_11 = lambda_x_[2*i],lambda_x_[2*i+1],lambda_y_[2*i],lambda_y_[2*i+1]
			a_00,a_01,a_10,a_11 = Tx_[3*i+1],Tx_[3*i+2],Ty_[3*i+1],Ty_[3*i+2]

			Tx_inv[4*i] = -a_11*tx + a_01*ty
			Tx_inv[4*i+1] = a_11 - lam_01*ty
			Tx_inv[4*i+2] = lam_11*tx - a_01
			Tx_inv[4*i+3] = - lam_11 + lam_01

			Ty_inv[4*i] = a_10*tx - a_00*ty
			Ty_inv[4*i+1] = lam_00*ty - a_10
			Ty_inv[4*i+2] = a_00 - lam_10*tx
			Ty_inv[4*i+3] = - lam_00 + lam_10

			lam_inv[4*i] = a_00*a_11 - a_10*a_01
			lam_inv[4*i+1] = a_10*lam_01 - a_11*lam_00
			lam_inv[4*i+2] = a_01*lam_10 - a_00*lam_11
			lam_inv[4*i+3] = lam_00*lam_11 - lam_01*lam_10

		Tx_inv[4] = Tx_inv[4] - (nx+WIDTH)*lam_inv[4]
		Tx_inv[5] = Tx_inv[5] - (nx+WIDTH)*lam_inv[5]
		Tx_inv[6] = Tx_inv[6] - (nx+WIDTH)*lam_inv[6]
		Tx_inv[7] = Tx_inv[7] - (nx+WIDTH)*lam_inv[7]

		Tx_inv[8] = Tx_inv[8] - (nx+WIDTH)*lam_inv[8]
		Tx_inv[9] = Tx_inv[9] - (nx+WIDTH)*lam_inv[9]
		Tx_inv[10] = Tx_inv[10] - (nx+WIDTH)*lam_inv[10]
		Tx_inv[11] = Tx_inv[11] - (nx+WIDTH)*lam_inv[11]

		Ty_inv[0] = Ty_inv[0] - (ny+WIDTH)*lam_inv[0]
		Ty_inv[1] = Ty_inv[1] - (ny+WIDTH)*lam_inv[1]
		Ty_inv[2] = Ty_inv[2] - (ny+WIDTH)*lam_inv[2]
		Ty_inv[3] = Ty_inv[3] - (ny+WIDTH)*lam_inv[3]

		Ty_inv[4] = Ty_inv[4] - (ny+WIDTH)*lam_inv[4]
		Ty_inv[5] = Ty_inv[5] - (ny+WIDTH)*lam_inv[5]
		Ty_inv[6] = Ty_inv[6] - (ny+WIDTH)*lam_inv[6]
		Ty_inv[7] = Ty_inv[7] - (ny+WIDTH)*lam_inv[7]

		transformations = numpy.row_stack((Tx_inv,Ty_inv,lam_inv))

		return transformations

	return None


def remap(directory, basename, tilts=None, mode=PROJECTIVE_TRANSFORMATIONS_1, qfix=[0,1], interpolation=NO_INTERPOLATION):

	filename = os.path.join(directory, basename + '_remap_stitch.mod')

	log.info("Remap model %s (with mode %i)." %(filename,mode))
	log.info("Interpolation: %i" %(interpolation))

	if not os.path.exists(filename):
		log.info('Model File %s for Stiching micrograhs does not exist!' %filename)
		raise SystemExit

	pts,nx,ny = readPoints(filename)

	frame = numpy.array([[1,2*nx],[1,2*ny]])

	resx = pts[:,0] - pts[:,3]
	resy = pts[:,1] - pts[:,4]

	res0 = numpy.mean(numpy.sqrt(resx**2+resy**2))

	if mode==TRANSLATIONS:	# Just translations

		(Tx,Ty,new_pts,res) = getPureTranslations(pts)

		log.info('Translations to apply:')
		for i in range(4):
			log.info('Quadrant %i' %(i+1))
			log.info('  %s' %Tx[3*i:3*i+3])
			log.info('  %s' %Ty[3*i:3*i+3])

		log.info('Original residual: %f' %res0)
		log.info('Residual after correction (min,max,mean): %s' %res)

	elif mode==AFFINE_TRANSFORMATIONS:	# for affine transformations

		(Tx,Ty,new_pts,res) = getAffineTransformations(pts,nx,ny)

		log.info('Affine Transformations to apply:')
		for i in range(4):
			log.info('Quadrant %i' %(i+1))
			log.info('  %s' %Tx[3*i:3*i+3])
			log.info('  %s' %Ty[3*i:3*i+3])

		log.info('Original residual: %f' %res0)
		log.info('Residual after correction (min,max,mean): %s' %res)

	elif mode==PROJECTIVE_TRANSFORMATIONS_1:	# for projective transformations

		(Tx,Ty,new_pts,res,frame) = getProjectiveTransformations_1(pts,nx,ny,qfix=qfix)

		log.info('Projective Transformations (v.1) to apply:')
		for i in range(4):
			log.info('Quadrant %i' %(i+1))
			log.info('  %s' %Tx[5*i:5*i+5])
			log.info('  %s' %Ty[5*i:5*i+5])

		log.info('Original residual: %f' %res0)
		log.info('Residual after correction (min,max,mean): %s' %res)

	elif mode==PROJECTIVE_TRANSFORMATIONS_2:	# for projective transformations

		(Tx,Ty,new_pts,res,frame) = getProjectiveTransformations_2(pts,nx,ny,qfix=qfix)

		log.info('Projective Transformations (v.2) to apply:')
		for i in range(4):
			log.info('Quadrant %i' %(i+1))
			log.info('  %s' %Tx[5*i:5*i+5])
			log.info('  %s' %Ty[5*i:5*i+5])

		log.info('Original residual: %f' %res0)
		log.info('Residual after correction (min,max,mean): %s' %res)

	# Calculate the transformtions. Inverse the map.

	transformations = inverseTransformation(Tx,Ty,nx,ny,mode)

	# Check the files

	if tilts==None or len(tilts)==0:
		tilts = []
		itilt = 0
		read_file = True
		while read_file:
			file_cam1 = os.path.join(directory, MRC, '%s_%03i_cam1.mrc' %(basename,itilt))
			file_cam2 = os.path.join(directory, MRC, '%s_%03i_cam2.mrc' %(basename,itilt))
			file_cam3 = os.path.join(directory, MRC, '%s_%03i_cam3.mrc' %(basename,itilt))
			file_cam4 = os.path.join(directory, MRC, '%s_%03i_cam4.mrc' %(basename,itilt))
			if os.path.exists(file_cam1) and os.path.exists(file_cam2) and \
			   os.path.exists(file_cam3) and os.path.exists(file_cam4):
				tilts.append(itilt)
			elif itilt==0:
				pass
			else:
				break
			itilt += 1

	ntilt = len(tilts)

	# Implement the remap

	xmin, xmax = frame[0].astype('int')
	ymin, ymax = frame[1].astype('int')

	final_file = os.path.join(directory,'%s_final_remap.mrc' %basename)

	f = mrc.MRCFile(final_file)
	f.setHeader(xmax-xmin+1,ymax-ymin+1,ntilt)

	final_section = numpy.zeros((2*nx,2*ny))

	index = 0

	for itilt in tilts:

		log.info('Tilt #%i' %itilt)

		t0 = time.time()

		final_section[:,:] = 0
		sections = []

		for i in range(1,5):
			file_ = os.path.join(directory, MRC, '%s_%03i_cam%i.mrc' %(basename,itilt,i))
			fq = mrc.MRCFile(file_)
			sections.append(fq.getZSliceAt(0))

		if mode==TRANSLATIONS or mode==AFFINE_TRANSFORMATIONS:
			misc.remap.remapSupp.remap_aff(sections,transformations,final_section)
		elif mode==PROJECTIVE_TRANSFORMATIONS_1 or mode==PROJECTIVE_TRANSFORMATIONS_2:
			misc.remap.remapSupp.remap_proj(sections, transformations, final_section, interpolation)

		t1 = time.time()

		f.setZSliceAt(index,final_section[xmin:xmax+1,ymin:ymax+1])

		index += 1

		t2 = time.time()

		log.info('Remap time: %.2f	Copy to MRC time: %.2f	Total Time: %.2f' %(t1-t0, t2-t1, t2-t0))


	f.updateHeader()

	os.system('imod %s' %f.filename)
