#  ===============================
#  GEORGIOS KARAGIANNIS (@PNNL.GOV)
#  GEORGIOS.KARAGIANNIS@PNNL.GOV
#  PNNL
#  2010-01-01
#  2013-08-20
#  ===============================
#
# GENERATE THE DATASET
#
# =============================================================================
# IMPORT BASIC MODULES 
# =============================================================================
#
import sys
from numpy import *

def vec(A) :
	return reshape(A,-1,order='F')

def system_evaluations(i) :
	#
	xi = array(loadtxt('x.dat')) ; 
	#
	u = array(loadtxt('y.dat')) ; u = array(u[:,i-1]) ;
	#
	x = array([1,2,3,4,5,6,7,8,9,10,11,12]) ; x = x[i-1] ;
	#
	return ( u, x, xi )

# =============================================================================
# MAIN PROGRAM
# =============================================================================

# GET xi,x, u =================================================================
#
id = 1
if len(sys.argv)>1 : id = int(sys.argv[1])
#
( u, x, xi ) = system_evaluations( id )
#
print( '------------------------') 
print(' GET xi,x, u ')
print( '------------------------') 
print( 'u.shape    : ' , u.shape ) 
print( 'u.mean     : ' , mean(u,0) ) 
print( 'u.std      : ' , std(u,0) ) 
print( 'x.shape    : ' , x.shape ) 
print( 'xi.shape   : ' , xi.shape) 
print( '------------------------\n') 
#
# SET THE PARAMETERS ==========================================================
#
xi_en = size(xi,0)                             # Size of the sample of y random input
xi_dim = size(xi,1)                             # Dimension of y random input
xi_min = xi.min()                             # 
xi_max = xi.max()                              # 
xi_degree = 4
#
print( '------------------------') 
print(' Algorithmic parameters : ')
print( '------------------------') 
print( 'x          : ' , x) 
print( 'xi_en      : ' , xi_en) 
print( 'xi_min     : ' , xi_min) 
print( 'xi_max     : ' , xi_max) 
print( 'xi_dim     : ' , xi_dim) 
print( 'xi_degree  : ' , xi_degree) 
print( '------------------------\n') 
#
# GET X, y ====================================================================
#
#from HermiteBases import DesignMatrixHermite
#from LegendreBases import DesignMatrixLegendre
from LegendreBases import DesignMatrixLegendre_transf
y = u
xi_min = array([ 0.01, 0.01, 0.1, 0.1, 9*10**(-7), 0.01, 0.01, 1e-6, 1e-6, 0.3 ])
xi_max = array([ 0.91, 1.0, 5.0, 5.0, 1e-1, 0.03, 21.0, 40000.0, 0.12, 0.6 ])

# xi_min = xi.min(0)
# xi_max = xi.max(0)
#( PSI, PSI_ind ) = DesignMatrixLegendre((2.0*xi-xi_max-xi_min)/(xi_max-xi_min), xi_en, xi_dim, xi_degree)
( PSI, PSI_ind ) = DesignMatrixLegendre_transf(xi, xi_en, xi_dim, xi_degree, xi_min, xi_max)
X = PSI
y_en = size(X,0)
dmax = size(X,1)
#
print( '------------------------') 
print(' GET X,y ')
print( '------------------------') 
print( 'y.shape    : ' , shape(y) ) 
print( 'X.shape    : ' , shape(X) ) 
print( 'y_en       : ' , y_en    ) 
print( 'dmax       : ' , dmax    ) 
print( '------------------------\n') 
#
# Save ========================================================================
#
print ( 'Saving...' )
#
import scipy.io as sio
#
# save raw data
#
save_file_name = './data/data.month='+str(id)+'.mat'
#
sio.savemat(\
			save_file_name, \
				{ \
					'x':x, \
					'xi':xi, \
					'xi_en':xi_en, \
					'xi_dim':xi_dim, \
					'xi_degree':xi_degree, \
					'xi_min':xi_min, \
					'xi_max':xi_max, \
					'y':y, \
					'X':X, \
					'y_en':y_en, \
					'X':X, \
					'dmax':dmax \
				} \
			)

