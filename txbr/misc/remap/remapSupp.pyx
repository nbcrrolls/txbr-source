cimport c_python as cpy
cimport c_numpy as cnp

cdef extern from "stdlib.h":
	ctypedef unsigned long size_t
	void *malloc(size_t size)
	void free(void *mem)

cnp.import_array()	# Numpy must be initialized

global NO_INTERPOLATION, BILINEAR_INTERPOLATION, BICUBIC_INTERPOLATION

NO_INTERPOLATION = 0
BILINEAR_INTERPOLATION = 1
BICUBIC_INTERPOLATION = 2

def remap_aff(sections, transformation, final_section):

	# Remap the four quadrants of one whole tilt

	cdef int nx, ny
	cdef double tx, axx, axy
	cdef double ty, ayx, ayy

	nx,ny = final_section.shape

	for iquadrant from 0<=iquadrant<4:

		tx, axx, axy = transformation[0,3*iquadrant:3*iquadrant+3]
		ty, ayx, ayy = transformation[1,3*iquadrant:3*iquadrant+3]

		remap_section_aff(iquadrant, sections[iquadrant], tx, axx, axy, ty, ayx, ayy, final_section, nx, ny)



def remap_proj(sections, transformation, final_section, interp):

	# Remap the four quadrants of one whole tilt

	cdef int nx, ny
	cdef double tx, a00, a01, a02
	cdef double ty, a10, a11, a12
	cdef double lam_0, lam_1, lam_01

	nx,ny = final_section.shape

	print transformation

	for iquadrant from 0<=iquadrant<4:

		print 'Quadrant #%i' %iquadrant
		print transformation[0,4*iquadrant:4*iquadrant+4]
		print transformation[1,4*iquadrant:4*iquadrant+4]
		print transformation[2,4*iquadrant:4*iquadrant+4]

		tx, a00, a01, a02 = transformation[0,4*iquadrant:4*iquadrant+4]
		ty, a10, a11, a12 = transformation[1,4*iquadrant:4*iquadrant+4]
		lam, lam_0, lam_1, lam_01 = transformation[2,4*iquadrant:4*iquadrant+4]

		remap_section_proj(iquadrant, sections[iquadrant], tx, a00, a01, a02, ty, a10, a11, a12, lam, lam_0, lam_1, lam_01, final_section, nx, ny, interp)



cdef remap_section_aff(int iquadrant, cnp.ndarray section, double tx, double axx, double axy, double ty, double ayx, double ayy,
		cnp.ndarray final_section, int nx, int ny):

	cdef int dbl_block

	dbl_block = 8	# 8 bytes to encode a double

	cdef int width, i
	cdef double *weight_x, *weight_y

	width = 1	# blend size outreach between the quadrants

	weight_x = <double *>malloc(dbl_block*nx)
	weight_y = <double *>malloc(dbl_block*ny)

	for i from 0<=i<nx: weight_x[i] = 0.0
	for i from 0<=i<ny: weight_y[i] = 0.0

	cdef int x0, x1, y0, y1
	cdef int nx_, ny_

	nx_ = nx/2
	ny_ = ny/2

	x0 = 1
	y0 = 1

	x1 = x0 + nx_
	y1 = y0 + ny_

	if iquadrant==0:

		x0 = 1
		x1 = x0 + nx_ + width
		for i from x0<=i<x1-2*width: weight_x[i] = 1.0
		for i from x1-2*width<=i<=x1: weight_x[i] = (x1-i)/(2.0*width)

		y0 = ny_ + 1 - width
		y1 = y0 + ny_
		for i from y0<=i<=y0+2*width: weight_y[i] = (i-y0)/(2.0*width)
		for i from y0+2*width<i<y1: weight_y[i] = 1.0

	elif iquadrant==1:

		x0 = nx_ + 1 - width
		x1 = x0 + nx_
		for i from x0<=i<=x0+2*width: weight_x[i] = (i-x0)/(2.0*width)
		for i from x0+2*width<i<x1: weight_x[i] = 1.0

		y0 = ny_ + 1 - width
		y1 = y0 + ny_
		for i from y0<=i<=y0+2*width: weight_y[i] = (i-y0)/(2.0*width)
		for i from y0+2*width<=i<y1: weight_y[i] = 1.0

	elif iquadrant==2:

		x0 = nx_ + 1 - width
		x1 = x0 + nx_
		for i from x0<=i<=x0+2*width: weight_x[i] = (i-x0)/(2.0*width)
		for i from x0+2*width<i<x1: weight_x[i] = 1.0

		y0 = 1
		y1 = y0 + ny_ + width
		for i from y0<=i<y1-2*width: weight_y[i] = 1.0
		for i from y1-2*width<=i<=y1: weight_y[i] = (y1-i)/(2.0*width)

	elif iquadrant==3:

		x0 = 1
		x1 = x0 + nx_ + width
		for i from x0<=i<=x1-2*width: weight_x[i] = 1.0
		for i from x1-2*width<i<=x1: weight_x[i] = (x1-i)/(2.0*width)

		y0 = 1
		y1 = y0 + ny_ + width
		for i from y0<=i<y1-2*width: weight_y[i] = 1.0
		for i from y1-2*width<=i<=y1: weight_y[i] = (y1-i)/(2.0*width)

	cdef int ix, iy, xv, yv
	cdef double value
	cdef double x_proj, y_proj, x_diff, y_diff

	cdef double *data, *data_q
	cdef int index, index_q_00, index_q_01, index_q_10, index_q_11
	cdef int strf_1, strf_2, strf_q_1, strf_q_2

	data  = <double *>cnp.PyArray_DATA(final_section)
	data_q  = <double *>cnp.PyArray_DATA(section)

	strf_1 = final_section.strides[0]/dbl_block
	strf_2 = final_section.strides[1]/dbl_block

	strf_q_1 = section.strides[0]/dbl_block
	strf_q_2 = section.strides[1]/dbl_block

	for ix from x0<=ix<x1:
		for iy from y0<=iy<y1:
			x_proj = tx + axx*ix + axy*iy
			y_proj = ty + ayx*ix + ayy*iy
			xv = <int>x_proj
			yv = <int>y_proj
			x_diff = x_proj - xv
			y_diff = y_proj - yv
			if xv>0 and xv<nx_-1 and yv>0 and yv<ny_-1:
				index = (ix-1)*strf_1 + (iy-1)*strf_2
				index_q_00 = xv*strf_q_1 + yv*strf_q_2
				index_q_01 = (xv)*strf_q_1 + (yv+1)*strf_q_2
				index_q_10 = (xv+1)*strf_q_1 + yv*strf_q_2
				index_q_11 = (xv+1)*strf_q_1 + (yv+1)*strf_q_2
				value = x_diff*y_diff*data_q[index_q_00] + (1.0-x_diff)*y_diff*data_q[index_q_10]	\
					+ x_diff*(1.0-y_diff)*data_q[index_q_01] + (1.0-x_diff)*(1.0-y_diff)*data_q[index_q_11]
				data[index] = data[index] + value*weight_x[ix]*weight_y[iy]

	free(weight_x)
	free(weight_y)



cdef remap_section_proj(int iquadrant, cnp.ndarray section, double tx, double a00, double a01, double a02,
		double ty, double a10, double a11, double a12, double lam, double lam_0, double lam_1, double lam_01,
		cnp.ndarray final_section, int nx, int ny, int interp):

	cdef int dbl_block

	dbl_block = 8	# 8 bytes to encode a double

	cdef int x0, x1, y0, y1
	cdef int nx_, ny_

	nx_ = nx/2
	ny_ = ny/2

	x0 = 1
	y0 = 1

	if iquadrant==0: x0, y0 = 1, ny_+1
	elif iquadrant==1: x0, y0 = nx_+1, ny_+1
	elif iquadrant==2: x0, y0 = nx_+1, 1
	elif iquadrant==3: x0, y0 = 1, 1

	x1 = x0 + nx_
	y1 = y0 + ny_

	cdef int ix, iy, xv, yv
	cdef double value
	cdef double x_proj, y_proj, x_diff, y_diff
	cdef double lambda_

	cdef double *data, *data_q
	cdef int index, index_q_00, index_q_01, index_q_10, index_q_11
	cdef int strf_1, strf_2, strf_q_1, strf_q_2

	data  = <double *>cnp.PyArray_DATA(final_section)
	data_q  = <double *>cnp.PyArray_DATA(section)

	strf_1 = final_section.strides[0]/dbl_block
	strf_2 = final_section.strides[1]/dbl_block

	strf_q_1 = section.strides[0]/dbl_block
	strf_q_2 = section.strides[1]/dbl_block

	cdef double secure_rate

	secure_rate = 0.1

	x0 = <int>((1-secure_rate)*x0)
	y0 = <int>((1-secure_rate)*y0)
	x1 = <int>((1+secure_rate)*x1)
	y1 = <int>((1+secure_rate)*y1)

	x0 = max(1,x0)
	y0 = max(1,y0)
	x1 = min(nx+1,x1)
	y1 = min(ny+1,y1)

	cdef int index_q_m1, index_q_0, index_q_1, index_q_2
	cdef double value_m1, value_0, value_1, value_2

	if interp==NO_INTERPOLATION:

		print 'x1=%i	y1=%i' %(x1,y1)

		for ix from x0<=ix<x1:

			for iy from y0<=iy<y1:

				x_proj = tx + a00*ix + a01*iy + a02*ix*iy
				y_proj = ty + a10*ix + a11*iy + a12*ix*iy
				lambda_ = lam + lam_0*ix + lam_1*iy + lam_01*ix*iy
				x_proj = x_proj/lambda_ - 1		# -1 for the C convention in the data array
				y_proj = y_proj/lambda_ - 1
				xv = <int>x_proj
				yv = <int>y_proj
				if xv>0 and xv<nx_-1 and yv>0 and yv<ny_-1:	# block
					index = (ix-1)*strf_1 + (iy-1)*strf_2
					index_q_0 = xv*strf_q_1 + yv*strf_q_2
					data[index] = data_q[index_q_0]


	if interp==BILINEAR_INTERPOLATION:

		for ix from x0<=ix<x1:

			for iy from y0<=iy<y1:

				x_proj = tx + a00*ix + (a01 + a02*ix)*iy
				y_proj = ty + a10*ix + (a11 + a12*ix)*iy
				lambda_ = lam + lam_0*ix + (lam_1 + lam_01*ix)*iy
				x_proj = x_proj/lambda_ - 1		# -1 for the C convention in the data array
				y_proj = y_proj/lambda_ - 1
				xv = <int>x_proj
				yv = <int>y_proj
				x_diff = x_proj - xv
				y_diff = y_proj - yv
				if xv>0 and xv<nx_-2 and yv>0 and yv<ny_-2:	# bilinear
					index = (ix-1)*strf_1 + (iy-1)*strf_2
					index_q_00 = xv*strf_q_1 + yv*strf_q_2
					index_q_01 = xv*strf_q_1 + (yv+1)*strf_q_2
					index_q_10 = (xv+1)*strf_q_1 + yv*strf_q_2
					index_q_11 = (xv+1)*strf_q_1 + (yv+1)*strf_q_2
					value = x_diff*y_diff*data_q[index_q_00] + x_diff*(1.0-y_diff)*data_q[index_q_01]	\
						  + (1.0-x_diff)*y_diff*data_q[index_q_10] + (1.0-x_diff)*(1.0-y_diff)*data_q[index_q_11]
					data[index] = value


	if interp==BICUBIC_INTERPOLATION:

		for ix from x0<=ix<x1:

			for iy from y0<=iy<y1:

				x_proj = tx + a00*ix + (a01 + a02*ix)*iy
				y_proj = ty + a10*ix + (a11 + a12*ix)*iy
				lambda_ = lam + lam_0*ix + (lam_1 + lam_01*ix)*iy
				x_proj = x_proj/lambda_ - 1		# -1 for the C convention in the data array
				y_proj = y_proj/lambda_ - 1
				xv = <int>x_proj
				yv = <int>y_proj
				x_diff = x_proj - xv
				y_diff = y_proj - yv

				if xv>1 and xv<nx_-3 and yv>1 and yv<ny_-3:	# cubic

					index = (ix-1)*strf_1 + (iy-1)*strf_2

					index_q_m1 = (xv-1)*strf_q_1 + (yv-1)*strf_q_2
					index_q_0 = xv*strf_q_1 + (yv-1)*strf_q_2
					index_q_1 = (xv+1)*strf_q_1 + (yv-1)*strf_q_2
					index_q_2 = (xv+2)*strf_q_1 + (yv-1)*strf_q_2
					value_m1 = p(x_diff, data_q[index_q_m1], data_q[index_q_0], data_q[index_q_1], data_q[index_q_2])

					index_q_m1 = (xv-1)*strf_q_1 + yv*strf_q_2
					index_q_0 = xv*strf_q_1 + yv*strf_q_2
					index_q_1 = (xv+1)*strf_q_1 + yv*strf_q_2
					index_q_2 = (xv+2)*strf_q_1 + yv*strf_q_2
					value_0 = p(x_diff, data_q[index_q_m1], data_q[index_q_0], data_q[index_q_1], data_q[index_q_2])

					index_q_m1 = (xv-1)*strf_q_1 + (yv+1)*strf_q_2
					index_q_0 = xv*strf_q_1 + (yv+1)*strf_q_2
					index_q_1 = (xv+1)*strf_q_1 + (yv+1)*strf_q_2
					index_q_2 = (xv+2)*strf_q_1 + (yv+1)*strf_q_2
					value_1 = p(x_diff, data_q[index_q_m1], data_q[index_q_0], data_q[index_q_1], data_q[index_q_2])

					index_q_m1 = (xv-1)*strf_q_1 + (yv+2)*strf_q_2
					index_q_0 = xv*strf_q_1 + (yv+2)*strf_q_2
					index_q_1 = (xv+1)*strf_q_1 + (yv+2)*strf_q_2
					index_q_2 = (xv+2)*strf_q_1 + (yv+2)*strf_q_2
					value_2 = p(x_diff, data_q[index_q_m1], data_q[index_q_0], data_q[index_q_1], data_q[index_q_2])

					value = p(y_diff, value_m1, value_0, value_1, value_2)
					data[index] = value



cdef p(double t, double a_m1, double a_0, double a_1, double a_2):

	# Centered differences
	return a_0 + t*(a_1 - a_m1 + t*(2.0*a_m1 - 5.0*a_0 + 4.0*a_1 - a_2 + t*(- a_m1 + 3.0*a_0 - 3.0*a_1 + a_2)))/2.0








