'''
Created on Feb 19, 2013

@author: georgios
CHECKED
'''

#  ===============================
#  GEORGIOS KARAGIANNIS (@PNNL.GOV)
#  GEORGIOS.KARAGIANNIS@PNNL.GOV
#  PNNL
#  2010-01-01
#  2013-08-20
#  ===============================
#

from numpy import *

# polymonials ==================================================================

def HermiteBases( xi, degree ) :
	#
	X = array( (degree+1)*[0.0] )
	X[0] = 1.0
	X[1] = xi 
	for j in range(2,degree+1) :
		X[j] = xi * X[j-1] -(j-1.0)*X[j-2]
	#
	return( X )

# Indeces =====================================================================

def gPC_bases_ind_rec( d , k ) : # recursive function
	#
	from scipy.misc import comb
	#
	i_row = 1
	#
	if d == 1 :
		A = array([[k]])
	else :
		xi = comb( k-1 , k-d , True )
		A = array( xi*[d*[0]] )
		for i in range(1,k-d+1+1) :
			xi = comb( k-i-1, k-i-d+1 , True )
			A[ i_row-1:i_row+xi-1 , 1-1]  = i
			A[ i_row-1:i_row+xi-1 , 2-1:d ] = gPC_bases_ind_rec(d-1, k-i)
			i_row = i_row + xi
	#
	return( A )

def gPC_bases_ind( d , k ) :
	A = gPC_bases_ind_rec(d,k+d) - 1
	return ( A )

def DesignMatrixHermite(xi, xi_en, xi_dim, degree, xi_mean, xi_std):
	#
	if degree == 0: 
		PSI = array(xi_en*[[1.0]])
		PSI_ind = array([[0]])
		return ( PSI, PSI_ind )
	#
	from scipy.misc import comb
	#
	dmax = comb(degree+xi_dim,xi_dim,True) 
	#
	# Compute the index matrix
	#
	PSI_ind = gPC_bases_ind( xi_dim , 0 ) 
	for k in range(1,degree+1) :
		PSI_ind = vstack( \
							[PSI_ind, 
							gPC_bases_ind( xi_dim , k )] \
						)
	#
	d = PSI_ind.shape[1]
	for j in range(1,int(d/2.0)+1):
		PSI_ind[:,[j-1,d+1-j-1]] = PSI_ind[:,[d+1-j-1,j-1]]
	#
	# Compute the design matrix
	#
	B = array( xi_dim*[xi_en*[(degree+1)*[0.0]]] )
	for j in range(1,xi_dim+1) :
		for i in range(1,xi_en+1) :
			zz = (xi[i-1,j-1]-xi_mean[j-1])/xi_std[j-1]
			B[j-1,i-1,0:degree+1] = HermiteBases(zz,degree).reshape(1,1,degree+1,order='F')
	#
	PSI = array( xi_en*[dmax*[0.0]] )
	for i in range(1,xi_en+1) :
		for j in range(1,dmax+1) :
			ind = copy(PSI_ind[j-1,0:xi_dim]) +1 
			c = 1
			for k in range( 1,xi_dim+1) :
				c = c*B[k-1,i-1,ind[k-1]-1] 
			PSI[i-1,j-1] = c 
	#
	return ( PSI, PSI_ind )

# Test me =====================================================================

if __name__ == "__main__":
	#
#    d = 9
#    k = 3
#    B = gPC_bases_ind( d , k )
#    print B
#    
#    xi = 2.1
#    degree = 3
#    X = LegendreBases( xi, degree )
#    print X
	#
	xi = array([[0.1, 0.2, 0.3, 0.4],[ 0.5, 0.6, 0.7, 0.8 ]])
	xi_en = 2
	xi_dim = 1
	degree = 3
	#
	[PSI, PSI_ind] = DesignMatrixHermite(xi, xi_en, xi_dim, degree)
	#
	print (PSI)
	print (PSI_ind)



