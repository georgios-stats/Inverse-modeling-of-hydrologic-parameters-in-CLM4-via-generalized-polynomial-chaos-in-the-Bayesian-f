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
def Energy(beta_mat, sig2_mat, uff, xi, xi_degree, xi_min, xi_max) :
    #
    xi_mu = array(1*[[0.561974,0.919636,2.813878,3.903571,0.000024,0.021146,19.453970,47.553009,0.003365,0.434089]])
    #
    #(PSI, PSI_ind) = DesignMatrixLegendre((2.0*xi-xi_max-xi_min)/(xi_max-xi_min), 1, size(xi,1), xi_degree)
    (PSI, PSI_ind) = DesignMatrixLegendre_transf(xi, 1, size(xi,1), xi_degree, xi_min, xi_max)
    mu = dot(PSI,beta_mat) 
    #
    En_Lik = 0.5*sum( ((uff-mu)**2)/sig2_mat ) 
    En_Pri = 0.5*sum( ((xi-xi_mu)**2) / 1000.0 ) \
                +sum( log(sig2_mat) )
    #
    logJ = sum( log(xi-xi_min) +log(xi_max-xi) -log(xi_max-xi_min) )
    #
    En = En_Lik +En_Pri -logJ
    #
    return En
#
def MH1(beta_mat, sig2_mat, uff, xi ,En, xi_degree, scl_hm1, xi_min, xi_max) :
    m = size(xi,1)
    #
    xi_N = xi.copy()
    xi_N = log(xi_N-xi_min) -log(xi_max-xi_N)
    xi_N = xi_N +scl_hm1*random.randn(m)
    xi_N = ( xi_min +xi_max*exp(xi_N) ) / (1.0+exp(xi_N))
    
    if sum(xi_N>=xi_max)>0 : return (xi, En, 0.0) 
    if sum(xi_N<=xi_min)>0 : return (xi, En, 0.0) 
    #
    En_N = Energy(beta_mat, sig2_mat, uff, xi_N, xi_degree, xi_min, xi_max) 
    acc = exp(min(0.0,-En_N+En)) 
    un = random.rand() 
    if (acc>=un) :
        xi_N = xi_N
        En = En_N
    else :
        xi_N = xi
        En_N = En
    #
    return (xi_N, En_N, acc) 
#
def MH2(beta_mat, sig2_mat, uff, xi ,En, xi_degree, scl_hm2, xi_min, xi_max) :
    #
    m = len(xi) ; 
    ee = scl_hm2*random.randn(m) 
    ee = ee/sqrt(sum(ee**2)) 
    #
    xi_N = xi.copy()
    xi_N = log(xi_N-xi_min) -log(xi_max-xi_N)
    xi_N = xi_N +scl_hm2*ee*random.randn()
    xi_N = ( xi_min +xi_max*exp(xi_N) ) / (1.0+exp(xi_N))
    #
    if sum(xi_N>=xi_max)>0 : return (xi, En, 0.0) 
    if sum(xi_N<=xi_min)>0 : return (xi, En, 0.0) 
    #
    En_N = Energy(beta_mat, sig2_mat, uff, xi_N, xi_degree, xi_min, xi_max) 
    acc = exp(min(0.0,-En_N+En)) 
    un = random.rand() 
    if (acc>=un) :
        xi_N = xi_N 
        En = En_N 
    else :
        xi_N = xi
        En_N = En
    #
    return (xi_N, En_N, acc) 
#
def updateSig2(a_sig2_mat,b_sig2_mat,beta_mat,uff,xi,xi_degree,xi_min,xi_max):
    m = len(uff)
    (PSI, PSI_ind) = DesignMatrixLegendre_transf(xi, 1, size(xi,1), xi_degree, xi_min, xi_max)
    mu = dot(PSI,beta_mat) 
    par1 = a_sig2_mat +0.5
    par2 = b_sig2_mat +0.5*(uff-mu)**2
    sig2_mat = random.gamma(par1,scale=1.0,size=(m,))
    sig2_mat = par2 / sig2_mat
    return sig2_mat
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
id = 1
#
N_sweep = 2*10**3
N_burnin = 10**3
scl_hm1 = 0.1 ;
scl_hm2 = 0.1 ;
#
a_sig2_mat = 10**(-3)
b_sig2_mat = 10**(-3)
#
file_name = './data/data.month='+str(id)+'.mat'
load_data = sio.loadmat(file_name)
xi_dim = int(load_data['xi_dim'])
xi_degree = int(load_data['xi_degree'])
xi_min = load_data['xi_min'].flatten() ;
xi_max = load_data['xi_max'].flatten() ;
#
file_name = './results/sample.bma.L2.month='+str(id)+'.mat'
load_data = sio.loadmat(file_name)
dmax = int(load_data['dmax'])
#
beta_mat = array( dmax*[uff_dim*[0.0]]) 
sig2_mat = array( 1*[uff_dim*[1.0]]) 
for id in range(1,uff_dim+1) :
    print ( 'Loading... ', id )
    file_name = './results/estimates.bma.L2.month='+str(id)+'.mat'
    load_data = sio.loadmat(file_name)
    beta_mat[:,id-1] = load_data['beta_est_bma_L2'].flatten() ;
    sig2_mat[:,id-1] = load_data['sig2_est_bma_L2'].flatten() ;
#    beta_mat[:,id-1] = load_data['beta_seed_bma_L2'].flatten() ;
#    sig2_mat[:,id-1] = load_data['sig2_seed_bma_L2'].flatten() ;
sig2_mat /= 10.0 ;
#
xi_chain = array( N_sweep*[xi_dim*[0.0]])
sig2_mat_chain = array( N_sweep*[uff_dim*[0.0]])
#xi = array(1*[xi_dim*[0.0]])
#xi = array([0.5*(xi_max+xi_min)])
xi = array(1*[[0.561974,0.919636,2.813878,3.903571,0.000024,0.021146,19.453970,47.553009,0.003365,0.434089]])
#
it_adapt_hm1 = 0 ; acc_adapt_hm1 = 0.0 ;
it_adapt_hm2 = 0 ; acc_adapt_hm2 = 0.0 ;
En = Energy(beta_mat, sig2_mat, uff, xi, xi_degree, xi_min, xi_max)
for it in range(-N_burnin,N_sweep) :
    #
    if False :
        (xi, En, acc_hm1) = MH1(beta_mat, sig2_mat, uff, xi ,En, xi_degree, scl_hm1, xi_min, xi_max) 
        it_adapt_hm1+=1
        acc_adapt_hm1+=acc_hm1
        gt = (1.0/it_adapt_hm1)**0.6
        scl_hm1 = exp( log(scl_hm1)+gt*(acc_hm1-0.234) )
    #
    if True :
        (xi, En, acc_hm2) = MH1(beta_mat, sig2_mat, uff, xi ,En, xi_degree, scl_hm2, xi_min, xi_max)
        it_adapt_hm2+=1 
        acc_adapt_hm2+=acc_hm2
        gt = (1.0/it_adapt_hm2)**0.6
        scl_hm2 = exp( log(scl_hm2)+gt*(acc_hm2-0.234) )
    #
    if False :
        sig2_mat = updateSig2(a_sig2_mat,b_sig2_mat,beta_mat,uff,xi,xi_degree,xi_min,xi_max)
    #
    if ( it>=0 ) : 
        xi_chain[it,:] = xi.flatten() ;
        sig2_mat_chain[it,:] = sig2_mat.flatten() ;
    #
#     print( it, acc_adapt_hm1/it_adapt_hm1, scl_hm1, acc_adapt_hm2/it_adapt_hm2, scl_hm2 )
#     print( it, acc_adapt_hm1/it_adapt_hm1, scl_hm1 )
#     print( it, acc_adapt_hm2/it_adapt_hm2, scl_hm2 )
#     if ( it>=0 ) :
#         print( mean(exp(xi_chain[0:it]),0), sqrt(var(exp(xi_chain[0:it]),0)/it) )
#    print ( (xi) )
#
if it_adapt_hm1>0 : acc_adapt_hm1 /= it_adapt_hm1
if it_adapt_hm2>0 : acc_adapt_hm2 /= it_adapt_hm2
#
xi_est = mean((xi_chain),0)
xi_se = sqrt( var((xi_chain),0)/N_sweep )
#
file_name = './results/sample.inv.bma.L2.month=all.mat'
sio.savemat( file_name, \
                { \
                 'sig2_mat_chain':sig2_mat_chain, \
                'xi_chain':xi_chain, \
                'xi_est':xi_est, \
                'xi_se':xi_se, \
                'acc_adapt_hm1':acc_adapt_hm1, \
                'scl_hm1':scl_hm1, \
                'acc_adapt_hm2':acc_adapt_hm2, \
                'scl_hm2':scl_hm2, \
                'uff':uff, \
                'uff_dim':uff_dim, \
                } \
            )
#
print(' ***END*** ')




