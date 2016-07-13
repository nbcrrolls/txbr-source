cimport c_python as cpy
cimport c_numpy as cnp

cdef extern from "stdlib.h":
	ctypedef int size_t
	ctypedef long intptr_t
	void *malloc(size_t size)
	void free(void* ptr)

cdef extern from "math.h":
	double sqrt(double)

cdef extern from "Python.h":
	ctypedef int Py_intptr_t

cdef struct PolyExp:
	int nnvar
	int nexp
	int *nterms
	int nterms_tot
	double *coeffs
	int *monoms


ONLY_FEATURE_DEPENDANT = 1
ONLY_TILT_DEPENDANT = 2
FEATURE_AND_TILT_DEPENDANT = 3
FEATURE_AND_TILT_INDEPENDANT = 4


import numpy
import logging
import time

cnp.import_array()	# Numpy must be initialized

cdef int nvar_core
cdef int *nvars_core
cdef int nnvar
cdef int npow_max_core
cdef int *npow_max_var_core
cdef PolyExp *E_core
cdef PolyExp *der_core
cdef PolyExp *hess_core
cdef int nder
cdef int nder_tot
cdef int *der_indices


#
# Evaluate the Error cost function, its derivatives and hessians
#
def initialize(core_parameters, ntilt, npatch):

	t0 = time.time()

	log = logging.getLogger("align")

	global nvar_core
	global nvars_core
	global nnvar
	global npow_max_core
	global npow_max_var_core
	global E_core
	global der_core
	global hess_core
	global nder
	global nder_tot
	global der_indices

	log.info('Global Calculations of Error, its derivatives and hessians')

	log.debug('Load core parameters')

	(nvars_core_,symbols,E_core_,nder_core_,modes_core_,der_core_,hess_core_) = core_parameters

	log.debug('Evaluate number and types of core variables:')

	cdef int i

	nvar_core = int(len(nvars_core_))

	nvars_core = <int*>malloc(sizeof(int)*nvar_core)

	nnvar = 0
	for i from 0<=i<nvar_core:
		nvars_core[i] = int(nvars_core_[i])
		log.debug('Number of Variable of Type %i: %i' %(i,nvars_core[i]))
		nnvar = nnvar + nvars_core[i]
		#nnvar += nvars_core[i]

	log.debug('Number of Core variable types: %i' %nvar_core)
	log.debug('Total Number of Core variables: %i' %nnvar)

	log.debug('Load Error Core expression parameters')

	E_core = initialize_core_Error(nnvar, E_core_)

	log.debug('Evaluate the maximum power each variable might experience in the Error expression!')

	npow_max_var_core = evaluate_max_power(nvar_core, nvars_core, E_core, &npow_max_core)

	log.debug('Maximum power for variable of any type: %i' %npow_max_core)

	log.debug('Power Time: %f' %(time.time()-t0))

	# Indexing -------------------------------------------------------------------------------------#

	cdef int ider

	nder = int(core_parameters[3])

	log.debug('Index the list of the core derivatives to all variables (any tilt/patch)...')

	der_indices = indexOfDerivatives(nder, modes_core_, ntilt, npatch, &nder_tot)

	log.debug('Total number of derivative variables: %i' %nder_tot)


	# Derivatives Matter ---------------------------------------------------------------------------#

	log.debug('Load the Derivatives Error Core expression parameters')

	tder = time.time()

	der_core = initialize_core_derivatives(nder,nnvar,der_core_)

	log.debug('%i different core derivatives expressions...' %der_core.nexp)

	log.debug('Total number of terms in the derivative expressions %i' %der_core.nterms_tot)

	log.debug('Derivatives Evaluation Timing: %f' %(time.time()-tder))

	# Hessians --------------------------------------------------------------------------------------#

	thess = time.time()

	log.debug('Load the Hessians Error Core expression parameters')

	hess_core = initialize_core_hessians(nder,nnvar,hess_core_)

	log.debug('%i different core hessians expressions...' %hess_core.nexp)

	log.debug('Total number of terms in the hessians expressions %i' %hess_core.nterms_tot)

	log.debug('Hessians Evaluation Timing: %f' %(time.time()-tder))

	log.debug('Initialization Time: %f' %(time.time()-t0))


#
# Evaluate the Error cost function, its derivatives and hessians
#
def values(core_parameters, ntilt, npatch, Emask, vars_):

	t0 = time.time()

	log = logging.getLogger('align')

	global nvar_core
	global nvars_core
	global nnvar
	global npow_max_core
	global npow_max_var_core
	global E_core
	global der_core
	global hess_core
	global nder
	global nder_tot
	global der_indices

	cdef int kmax
	kmax = npow_max_core + 1

	log.debug('Convert all variables')

	#print 'nvar_core=%i' %nvar_core
	#print 'l=%i' %len(vars_)

	if nvar_core!=len(vars_):
		raise ValueError('Numbers of variable types do not match!')

	cdef int nnvar_tot

	nnvar_tot = 0
	_vars_ = []

	for ivar from 0<=ivar<nvar_core:
		l = len(vars_[ivar][0])
		log.debug('Number of type %i variables: %i' %(ivar,l))
		nnvar_tot = nnvar_tot + l
		_vars_ = _vars_ + vars_[ivar][0]

	log.debug('Total Number of variables %i' %nnvar_tot)

	cdef int *modes

	modes = <int*>malloc(sizeof(int)*nvar_core)

	for ivar from 0<=ivar<nvar_core:
		modes[ivar] = int(vars_[ivar][1])

	cdef int *skips

	skips = <int*>malloc(sizeof(int)*nnvar_tot)

	cdef int index
	index = 0

	for ivar from 0<=ivar<nvar_core:
		#sk =  vars_[ivar][2]
		#l = len(sk)
		l = len(vars_[ivar][0])
		for i from 0<=i<l:
			#skips[index] = int(sk[i])
			skips[index] = 0
			index = index + 1

	log.debug('Prepare a specific table of variables for the variables')

	cdef double *vars_pow

	vars_pow = prepare_variables(nvar_core, nvars_core, modes, skips, npow_max_core, npow_max_var_core, ntilt, npatch, nnvar_tot, _vars_)

	log.debug('Call the C Routines')

	t1 = time.time()

	Err = error(nnvar, E_core, ntilt, npatch, Emask, kmax, vars_pow)

	t2 = time.time()

	Der = derivatives(nder, nnvar, der_core, ntilt, npatch, Emask, skips, nder_tot, der_indices, kmax, vars_pow)

	t3 = time.time()

	Hess = hessians(nder, nnvar, hess_core, ntilt, npatch, Emask, skips, nder_tot, der_indices, kmax, vars_pow)

	t4 = time.time()

	log.debug('Timing: %f	Err %f	Der %f	Hess %f' %(t4-t1,t2-t1,t3-t2,t4-t3))

	free(modes)
	free(skips)
	free(vars_pow)

	log.debug('Total Evaluation Time=%f' %(time.time()-t0))

	return (Err,Der,Hess)


#
# Evaluate the maximum power that each variable type will experience in the Error cost function
#
cdef int *evaluate_max_power(int nvar_core, int *nvars_core, PolyExp *E_core, int *npow_max_core):

	cdef int nterms_core
	cdef int *monoms_core

	nterms_core = E_core.nterms_tot
	monoms_core = E_core.monoms

	# npow_max_core: Will be assigned the maximum power that any variable might experience

	cdef int* npow_max_var_core
	cdef int ivar,index

	npow_max_var_core = <int*>malloc(sizeof(int)*nvar_core)

	index = 0
	for ivar from 0<=ivar<nvar_core:
		npow_max_var_core[ivar] = 0
		for ii from 0<=ii<nvars_core[ivar]:
			for k from 0<=k<nterms_core:
				npow_max_var_core[ivar] = max(npow_max_var_core[ivar],monoms_core[index])
				index = index + 1

	npow_max_core[0] = 0	# For every variable type

	for ivar from 0<=ivar<nvar_core:
		if npow_max_core[0]<npow_max_var_core[ivar]:
			npow_max_core[0] = npow_max_var_core[ivar]

	return npow_max_var_core

#
# Prepare the variables before calculating the Error cost function and its derivatives/hessians
#
cdef double *prepare_variables(int nvar_core, int *nvars_core, int *modes, int *skips, int npow_max_core, int* npow_max_var_core, int ntilt, int npatch, int nnvar_tot, object vars):

	# nvar_core: Number of core variable types (usually 4): a, b, c, d ....
	# nvars_core: Array of nvar int, containing the number of core variables for each type
	# npow_max: Maximum power at which a variable will be used
	# npow_max_var_core: Array of int indicating the maximum power that a variable of a given type will be at.
	# ntilt: Number of tilts
	# npatch: Number of patches
	# vars: Array of the variables

	cdef int itilt,ipatch,index,ivar,i0,ii,k,nnvar,nmax,i_,i
	cdef double *vars_pow
	cdef int *offset

	# Evaluate the offsets for each type of variables

	offset = <int*>malloc(sizeof(int)*nvar_core)

	offset[0] = 0
	for ivar from 1<=ivar<nvar_core:
		offset[ivar] = offset[ivar-1] + nvars_core[ivar-1]

	# Allocate memory for the buffer containing the powered items

	nnvar = 0
	for ivar from 0<=ivar<nvar_core:
		nnvar = nnvar + nvars_core[ivar]

	# Build the index array

	cdef int *indices

	indices = <int*>malloc(sizeof(int)*ntilt*npatch*nnvar)

	index = 0
	i_ = 0

	for ivar from 0<=ivar<nvar_core:
		for i from 0<=i<nvars_core[ivar]:
			if modes[ivar]==ONLY_FEATURE_DEPENDANT:
				for ipatch from 0<=ipatch<npatch:
					for itilt from 0<=itilt<ntilt:
						indices[itilt*npatch*nnvar+ipatch*nnvar+i_] = index
					index = index + 1
			if modes[ivar]==ONLY_TILT_DEPENDANT:
				for itilt from 0<=itilt<ntilt:
					for ipatch from 0<=ipatch<npatch:
						indices[itilt*npatch*nnvar+ipatch*nnvar+i_] = index
					index = index + 1
			if modes[ivar]==FEATURE_AND_TILT_DEPENDANT:
				for itilt from 0<=itilt<ntilt:
					for ipatch from 0<=ipatch<npatch:
						indices[itilt*npatch*nnvar+ipatch*nnvar+i_] = index
						index = index + 1
			if modes[ivar]==FEATURE_AND_TILT_INDEPENDANT:
				for itilt from 0<=itilt<ntilt:
					for ipatch from 0<=ipatch<npatch:
						indices[itilt*npatch*nnvar+ipatch*nnvar+i_] = index
				index = index + 1
			i_ = i_ + 1

#		print "ivar: %i->%i  n: %i  mode: %i" %(ivar,index,nvars_core[ivar],modes[ivar])


	# Calculate the powered variables

	vars_pow = <double*>malloc(sizeof(double)*(npow_max_core+1)*nnvar*ntilt*npatch)

	for itilt from 0<=itilt<ntilt:

		for ipatch from 0<=ipatch<npatch:

			i_ = 0

			for ivar from 0<=ivar<nvar_core:

				nmax = npow_max_var_core[ivar]+1

				for ii from 0<=ii<nvars_core[ivar]:

					index = indices[(itilt*npatch+ipatch)*nnvar+i_]


					i0 = (itilt*npatch+ipatch)*(npow_max_core+1)*nnvar + (offset[ivar]+ii)*(npow_max_core+1)

					#print 'index=%i/%i	itilt=%i  ipatch=%i  i0=%i   max=%i' %(index,len(vars),itilt,ipatch,i0,(npow_max_core+1)*nnvar*ntilt*npatch)

					#if skips[i_]==1:
					#if skips[index]==1:
					if skips[index]==-1:	# never happen

						vars_pow[i0] = 1.0
						for k from 1<=k<nmax: vars_pow[i0+k] = 0.0

					else:

						vars_pow[i0] = 1.0
						vars_pow[i0+1] = vars[index]

						for k from 2<=k<nmax: vars_pow[i0+k] = vars_pow[i0+(k-1)]*vars_pow[i0+1]

				#	if ipatch==0 and itilt==0:
				#		print 'ivar=%i ii=%i i0=%i index=%i  var=%f' %(ivar,ii,i0,index,vars[index])

					i_ = i_ + 1


	free(offset)
	free(indices)

#	print vars[0]


#	index = 0
#	i_=0
#	kmax = npow_max_core+1
#	for ivar from 0<=ivar<nvar_core:
#		nmax = npow_max_var_core[ivar]+1
#		for ii from 0<=ii<nvars_core[ivar]:
#			for k from 0<=k<nmax:
#				print 'ivar=%-5i ii=%-5i i=%-5i k=%-5i index=%-5i var=%10.4e	var_pow=%10.4e' %(ivar,ii,index,k,index,vars[index],vars_pow[index*kmax+k])
#			index = index + 1

	return vars_pow

#
# Map the derivative indices
#
cdef int *indexOfDerivatives(int nder, object modes, int ntilt, int npatch, int *nder_tot):

	cdef int itilt,ipatch,index

	nder_tot[0] = 0

	for ider from 0<=ider<nder:
		if modes[ider]==ONLY_FEATURE_DEPENDANT: nder_tot[0] = nder_tot[0] + npatch
		if modes[ider]==ONLY_TILT_DEPENDANT: nder_tot[0] = nder_tot[0] + ntilt
		if modes[ider]==FEATURE_AND_TILT_DEPENDANT: nder_tot[0] = nder_tot[0] + ntilt*npatch
		if modes[ider]==FEATURE_AND_TILT_INDEPENDANT: nder_tot[0] = nder_tot[0] + 1

	cdef int *der_indices

	der_indices = <int*>malloc(sizeof(int)*ntilt*npatch*nder)

	index = 0

	for ider from 0<=ider<nder:
		if modes[ider]==ONLY_FEATURE_DEPENDANT:
			for ipatch from 0<=ipatch<npatch:
				for itilt from 0<=itilt<ntilt:
					der_indices[itilt*npatch*nder+ipatch*nder+ider] = index
				index = index + 1
		if modes[ider]==ONLY_TILT_DEPENDANT:
			for itilt from 0<=itilt<ntilt:
				for ipatch from 0<=ipatch<npatch:
					der_indices[itilt*npatch*nder+ipatch*nder+ider] = index
				index = index + 1
		if modes[ider]==FEATURE_AND_TILT_DEPENDANT:
			for itilt from 0<=itilt<ntilt:
				for ipatch from 0<=ipatch<npatch:
					der_indices[itilt*npatch*nder+ipatch*nder+ider] = index
					index = index + 1
		if modes[ider]==FEATURE_AND_TILT_INDEPENDANT:
			for itilt from 0<=itilt<ntilt:
				for ipatch from 0<=ipatch<npatch:
					der_indices[itilt*npatch*nder+ipatch*nder+ider] = index
			index = index + 1

	return der_indices


cdef PolyExp *initialize_core_Error(int nnvar, object E_core_):

	(E_coeffs_core_, E_monoms_core_) = E_core_

	cdef int *E_nterms_core,*E_monoms_core
	cdef double *E_coeffs_core
	cdef int k

	E_nterms_core = <int*>malloc(sizeof(int))

	E_nterms_core[0] = int(len(E_coeffs_core_))

	E_coeffs_core = <double*>malloc(sizeof(double)*E_nterms_core[0])

	for k in range(E_nterms_core[0]):
		E_coeffs_core[k] = float(E_coeffs_core_[k])

	E_monoms_core = <int*>malloc(sizeof(int)*E_nterms_core[0]*nnvar)

	cdef index

	index = 0

	for i from 0<=i<nnvar:
		for k from 0<=k<E_nterms_core[0]:
			E_monoms_core[index] = int(E_monoms_core_[k][i])
			index = index + 1

	cdef PolyExp *E_core

	E_core = <PolyExp*>malloc(sizeof(PolyExp))

	E_core.nnvar = nnvar
	E_core.nexp = 1
	E_core.nterms = E_nterms_core
	E_core.nterms_tot = E_nterms_core[0]
	E_core.coeffs = E_coeffs_core
	E_core.monoms = E_monoms_core

	return E_core



cdef PolyExp *initialize_core_derivatives(int nder, int nnvar, object der_core_):

	cdef int der_nterms_core_tot, *der_nterms_core, *der_monoms_core
	cdef cnp.ndarray der_coeffs_core_,der_monoms_core_
	cdef double *der_coeffs_core

	der_nterms_core_tot = 0

	der_nterms_core = <int*>malloc(sizeof(int)*nder)

	for ider from 0<=ider<nder:
		der_nterms_core[ider] = int(len(der_core_[ider][0]))
		der_nterms_core_tot = der_nterms_core_tot + der_nterms_core[ider]

	der_coeffs_core = <double*>malloc(sizeof(double)*der_nterms_core_tot)
	der_monoms_core = <int*>malloc(sizeof(int)*der_nterms_core_tot*nnvar)

	index = 0

	for ider from 0<=ider<nder:
		der_coeffs_core_ = der_core_[ider][0]
		for k in range(der_nterms_core[ider]):
			der_coeffs_core[index] = float(der_coeffs_core_[k])
			index = index + 1

	index = 0

	for ider from 0<=ider<nder:
		der_monoms_core_ = der_core_[ider][1]
		for i from 0<=i<nnvar:
			for k from 0<=k<der_nterms_core[ider]:
				der_monoms_core[index] = int(der_monoms_core_[k][i])
				index = index + 1

	cdef PolyExp *der_core

	der_core = <PolyExp*>malloc(sizeof(PolyExp))

	der_core.nnvar = nnvar
	der_core.nexp = nder
	der_core.nterms_tot = der_nterms_core_tot
	der_core.nterms = der_nterms_core
	der_core.coeffs = der_coeffs_core
	der_core.monoms = der_monoms_core

	return der_core



cdef PolyExp *initialize_core_hessians(int nder, int nnvar, object hess_core_):

	cdef int ihess0,ihess1,ihess2,nhess_,nhess
	cdef int hess_nterms_core_tot, *hess_nterms_core, *hess_monoms_core
	cdef double *hess_coeffs_core

	#nhess_ = int(len(hess_core_))
	nhess_ = nder
	nhess = int(nhess_*(nhess_+1.0)/2.0)

	hess_nterms_core_tot = 0

	hess_nterms_core = <int*>malloc(sizeof(double)*nhess)

	index = 0

	for ihess1 from 0<=ihess1<nhess_:
		for ihess2 from ihess1<=ihess2<nhess_:
			hess_coeffs_core_ = hess_core_[ihess1][ihess2][0]
			hess_nterms_core[index] = int(len(hess_coeffs_core_))
			hess_nterms_core_tot = hess_nterms_core_tot + hess_nterms_core[index]
			index = index + 1


	hess_coeffs_core = <double*>malloc(sizeof(double)*hess_nterms_core_tot)

	index = 0
	for ihess1 from 0<=ihess1<nhess_:
		for ihess2 from ihess1<=ihess2<nhess_:
			hess_coeffs_core_ = hess_core_[ihess1][ihess2][0]
			ihess0 = (ihess1*(2*nhess_-ihess1-1) + 2*ihess2)/2
			for k in range(hess_nterms_core[ihess0]):
				hess_coeffs_core[index] = float(hess_coeffs_core_[k])
				index = index + 1

	tt = time.time()

	hess_monoms_core = <int*>malloc(sizeof(int)*hess_nterms_core_tot*nnvar)

	index = 0
	for ihess1 from 0<=ihess1<nhess_:
		for ihess2 from ihess1<=ihess2<nhess_:
			hess_monoms_core_ = hess_core_[ihess1][ihess2][1]
			ihess0 = (ihess1*(2*nhess_-ihess1-1) + 2*ihess2)/2
			for i from 0<=i<nnvar:
				for k from 0<=k<hess_nterms_core[ihess0]:
					hess_monoms_core[index] = int(hess_monoms_core_[k][i])
					index = index + 1

	cdef PolyExp *hess_core

	hess_core = <PolyExp*>malloc(sizeof(PolyExp))

	hess_core.nnvar = nnvar
	hess_core.nexp = nhess
	hess_core.nterms_tot = hess_nterms_core_tot
	hess_core.nterms = hess_nterms_core
	hess_core.coeffs = hess_coeffs_core
	hess_core.monoms = hess_monoms_core

	return hess_core



#
# Calculate the Error cost function. The variables should have been previously ordered
# in a specific array in the prepare_variables(...) subroutine
#
cdef double error(int nnvar, PolyExp *E_core, int ntilt, int npatch, cnp.ndarray Emask, int kmax, double *vars_pow):

	# nnvar: Number of variables (of all types) in the Error core expression
	# nterms_core: Number of coefficients and monoms in the Error core expression
	# coeffs_core: Array of nterms_core double for the Error core expression coefficients
	# monoms_core: Array of nterms_core*nnvar double for the Error core expression monoms
	# ntilt: Number of tilts
	# npatch: Number of patches
	# kmax: Maximum number of powered items for any variable
	# vars_pow: Array of the powered variables (see prepare_variables functions)

	cdef int nterms_core
	cdef double* coeffs_core
	cdef int* monoms_core

	nterms_core = E_core.nterms_tot
	coeffs_core = E_core.coeffs
	monoms_core = E_core.monoms

	cdef int itilt,ipatch,k,ii,index
	cdef double Err,*Err_

	cdef char* data_mask
	cdef double* mask

	data_mask = Emask.data

	Err = 0

	Err_ = <double*>malloc(sizeof(double)*nterms_core)	# will be initialized in the loop

	cdef e

	index = 0
	for itilt from 0<=itilt<ntilt:
		for ipatch from 0<=ipatch<npatch:
			mask = <double*>(data_mask + itilt*Emask.strides[0] + ipatch*Emask.strides[1])
			for k from 0<=k<nterms_core:
				Err_[k] = mask[0]*coeffs_core[k]
			for ii from 0<=ii<nnvar:
				for k from 0<=k<nterms_core:
					Err_[k] = Err_[k]*vars_pow[index+monoms_core[nterms_core*ii+k]]
				#	if vars_pow[index+monoms_core[nterms_core*ii+k]]!=1:
				#		print 'ii=%i/%i k=%i/%i  coeff=%i  monom=%i  var=%f' %(ii,nnvar,k,nterms_core,coeffs_core[k],monoms_core[nterms_core*ii+k],vars_pow[index+monoms_core[nterms_core*ii+k]])
				#if itilt==1 and ipatch==0: print 'coeff in pyrex: %f' %vars_pow[index+1]
				#if ii>=nnvar-2: print 'coeff in pyrex: %f' %vars_pow[index+1]
				index = index + kmax
			e = 0
			for k from 0<=k<nterms_core:
				Err = Err + Err_[k]
				e = e + Err_[k]
	#		print 'itilt=%i ipatch=%i  Err=%f' %(itilt,ipatch,e)

	#print 'total error=%f' %Err

	free(Err_)

	return Err


#
# Calculate the derivatives of the cost function
#
cdef cnp.ndarray derivatives(int nder, int nnvar, PolyExp *der_core, int ntilt, int npatch, cnp.ndarray Emask, int *skips, int nder_tot, int *der_indices, int kmax, double *vars_pow):

	# nder: Total number of core derivatives
	# nnvar: Total number of variables in the derivatives core expression (same as in the Error core function)
	# *nterms_core: Array containing the number of coefficients and monoms in the derivatives core expressions
	# coeffs_core: Array of nder*nterms_core double for the derivatives expression coefficients
	# monoms_core: Array of nder*nterms_core*nnvar double for the derivatives expression monoms
	# ntilt: Number of tilts
	# npatch: Number of patches
	# kmax: Maximum number of powered items for any variable
	# vars_pow: Array of the powered variables (see prepare_variables functions)

	cdef int *nterms_core
	cdef double *coeffs_core
	cdef int *monoms_core

	nterms_core = der_core.nterms
	coeffs_core = der_core.coeffs
	monoms_core = der_core.monoms

	cdef int ider,itilt,ipatch,icoeff,k,ii,index,index0
	cdef nterms_core_max
	cdef double *der_

	nterms_core_max = 0
	for ider from 0<=ider<nder:
		if nterms_core_max<nterms_core[ider]:
			nterms_core_max = nterms_core[ider]

	cdef char* data_mask
	cdef double* mask

	data_mask = Emask.data

	cdef cnp.ndarray Der
	cdef double *value
	cdef char* data

	Der = numpy.zeros((nder_tot),dtype=float)
	data = Der.data

	der_ = <double*>malloc(sizeof(double)*nterms_core_max)	# will be initialized in the loop

	icoeff = 0

	for ider from 0<=ider<nder:
		index = 0
		for itilt from 0<=itilt<ntilt:
			for ipatch from 0<=ipatch<npatch:
				mask = <double*>(data_mask + itilt*Emask.strides[0] + ipatch*Emask.strides[1])
				index0 = der_indices[itilt*npatch*nder+ipatch*nder+ider]
				value = <double*>(data + index0*Der.strides[0])
				for k from 0<=k<nterms_core[ider]:
					der_[k] = mask[0]*coeffs_core[icoeff+k]
				for ii from 0<=ii<nnvar:
					for k from 0<=k<nterms_core[ider]:
						der_[k] = der_[k] *vars_pow[index+monoms_core[icoeff*nnvar+nterms_core[ider]*ii+k]]
					index = index + kmax
				for k from 0<=k<nterms_core[ider]:
					value[0] = value[0] + <double>der_[k]
		icoeff = icoeff + nterms_core[ider]

	free(der_)

	for ider from 0<=ider<nder_tot:
		if skips[ider]==1: Der[ider] = 0.0

	return Der

#
# Calculate the hessians of the cost function
#
cdef cnp.ndarray hessians(int nder, int nnvar, PolyExp *hess_core, int ntilt, int npatch, cnp.ndarray Emask, int *skips, int nder_tot, int *der_indices, int kmax, double *vars_pow):

	# nhess: Number of core derivation variables involved in hessians. Number of hessians expressions: nhess*nhess
	# nnvar: Total number of variables in the hessians core expression (same as in the Error core function)
	# *nterms_core: Array containing the number of coefficients and monoms in the hessians core expressions
	# coeffs_core: Array of nhess*nterms_core double for the hessians expression coefficients
	# monoms_core: Array of nhess*nterms_core*nnvar double for the hessians expression monoms
	# ntilt: Number of tilts
	# npatch: Number of patches
	# kmax: Maximum number of powered items for any variable
	# vars_pow: Array of the powered variables (see prepare_variables functions)


	cdef int *nterms_core
	cdef double *coeffs_core
	cdef int *monoms_core

	nterms_core = hess_core.nterms
	coeffs_core = hess_core.coeffs
	monoms_core = hess_core.monoms

	cdef int nhess_,nhess,nterms_core_max
	cdef int ihess,ihess1,ihess2,itilt,ipatch,icoeff,k,ii,index,index1,index2
	cdef double *hess_

	nhess_ = nder
	nhess = int(nder*(nder+1)/2)

	#nhess_ = int((-1.0+sqrt(1.0+8.0*nhess))/2.0)

	nterms_core_max = 0
	for ihess from 0<=ihess<nhess:
		if nterms_core_max<nterms_core[ihess]:
			nterms_core_max = nterms_core[ihess]

	cdef cnp.ndarray Hess
	cdef double *value1,*value2
	cdef char* data

	cdef char* data_mask
	cdef double* mask

	data_mask = Emask.data

	Hess = numpy.zeros((nder_tot,nder_tot),dtype=float)
	data = Hess.data

	hess_ = <double*>malloc(sizeof(double)*nterms_core_max)	# will be initialized in the loop

	icoeff = 0

	for ihess1 from 0<=ihess1<nhess_:
		for ihess2 from ihess1<=ihess2<nhess_:
			ihess = (ihess1*(2*nhess_-ihess1-1)+2*ihess2)/2
			index = 0
			for itilt from 0<=itilt<ntilt:
				for ipatch from 0<=ipatch<npatch:
					mask = <double*>(data_mask + itilt*Emask.strides[0] + ipatch*Emask.strides[1])
					index1 = der_indices[itilt*npatch*nhess_+ipatch*nhess_+ihess1]
					index2 = der_indices[itilt*npatch*nhess_+ipatch*nhess_+ihess2]
					value1 = <double*>(data + index1*Hess.strides[0] + index2*Hess.strides[1])
					value2 = <double*>(data + index2*Hess.strides[0] + index1*Hess.strides[1])
					for k from 0<=k<nterms_core[ihess]:
						hess_[k] = mask[0]*coeffs_core[icoeff+k]
					for ii from 0<=ii<nnvar:
						for k from 0<=k<nterms_core[ihess]:
							hess_[k] = hess_[k] *vars_pow[index+monoms_core[icoeff*nnvar+nterms_core[ihess]*ii+k]]
						index = index + kmax
					for k from 0<=k<nterms_core[ihess]:
						value1[0] = value1[0] + <double>hess_[k]
						if index1!=index2: # Take the symetric contribution
							value2[0] = value2[0] + <double>hess_[k]
			icoeff = icoeff +nterms_core[ihess]

	free(hess_)

	for ihess1 from 0<=ihess1<nder_tot:
		for ihess2 from 0<=ihess2<nder_tot:
			if skips[ihess1]==1: Hess[ihess1,ihess2] = 0.0
			if skips[ihess2]==1: Hess[ihess1,ihess2] = 0.0

	return Hess


