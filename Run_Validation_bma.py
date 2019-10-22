#  ===============================
#  GEORGIOS KARAGIANNIS (@PNNL.GOV)
#  GEORGIOS.KARAGIANNIS@PNNL.GOV
#  PNNL
#  2013-08-20
#  2014-09-15
#  ===============================
#
# =============================================================================
# IMPORT BASIC MODULES 
# =============================================================================
#
import sys
from numpy import *
#from LegendreBases import DesignMatrixLegendre
from LegendreBases import DesignMatrixLegendre_transf
#
def vec(A) : 
    return reshape(A,-1,order='F')
#
def gPCE(beta_mat,xi, xi_degree, xi_min, xi_max) :
    #
    #(PSI, PSI_ind) = DesignMatrixLegendre((2.0*xi-xi_max-xi_min)/(xi_max-xi_min), 1, size(xi,1), xi_degree)
    (PSI, PSI_ind) = DesignMatrixLegendre_transf(xi, 1, size(xi,1), xi_degree, xi_min, xi_max)
    mu = dot(PSI,beta_mat) 
    #
    return mu
#
# =============================================================================
# MAIN PROGRAM
# =============================================================================
#
# GET beta_mat, sig2_mat ======================================================
#
import scipy.io as sio
#
#uff = array([[5.208, 8.837, 21.771, 43.011, 56.944, 69.185, 70.860, 51.585, 39.797, 16.491, 9.046, 5.245]]) ;
uff = array([[15.542,22.017,41.365,59.095,58.377,58.813,45.107,41.362,31.250,28.645,17.635,12.778]]) ;
uff_dim = size(uff,1)
#
idd = 1
#
N_sweep = 2*10**3
#
file_name = './data/data.month='+str(idd)+'.mat'
load_data = sio.loadmat(file_name)
xi_dim = int(load_data['xi_dim'])
xi_degree = int(load_data['xi_degree'])
xi_min = load_data['xi_min'].flatten() ;
xi_max = load_data['xi_max'].flatten() ;
#
file_name = './results/sample.bma.L2.month='+str(idd)+'.mat'
load_data = sio.loadmat(file_name)
dmax = int(load_data['dmax'])
#
beta_mat = array( dmax*[uff_dim*[0.0]]) 
sig2_mat = array( 1*[uff_dim*[0.0]]) 
for idd in range(1,uff_dim+1) :
    print ( 'Loading... ', idd )
    file_name = './results/estimates.bma.L2.month='+str(idd)+'.mat'
    load_data = sio.loadmat(file_name)
    beta_mat[:,idd-1] = load_data['beta_est_bma_L2'].flatten() ;
    sig2_mat[:,idd-1] = load_data['sig2_est_bma_L2'].flatten() ;
#    beta_mat[:,idd-1] = load_data['beta_seed_bma_L2'].flatten() ;
#    sig2_mat[:,idd-1] = load_data['sig2_seed_bma_L2'].flatten() ;
#
file_name = './results/sample.inv.bma.L2.month=all.mat'
load_data = sio.loadmat(file_name)
xi_chain = array(load_data['xi_chain']) ;
#
uff_chain = array( N_sweep*[uff_dim*[0.0]])
#
for t in range(1,N_sweep+1) :
    xi = array([xi_chain[t-1,:].copy()])
    uff_chain[t-1,:] = gPCE(beta_mat,xi, xi_degree, xi_min, xi_max)
#
file_name = './results/sample.valid.bma.L2.month=all.mat'
sio.savemat( file_name, \
                { \
                'uff_chain':uff_chain, \
                'xi_chain':xi_chain, \
                'N_sweep':N_sweep, \
                'xi_degree':xi_degree, \
                'xi_min':xi_min, \
                'xi_max':xi_max \
                } \
            )
#
print(' ***END*** ')






