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
#
def vec(A) : 
	return reshape(A,-1,order='F')
#
# =============================================================================
# MAIN PROGRAM
# =============================================================================
#
id = 1
if len(sys.argv)>1 : id = int(sys.argv[1])
#
# LOAD DATA ===================================================================
#
print ( 'Loading... ', id )
#
import scipy.io as sio
#
# load_file_name = './data/data.'+str(x_en)+'.'+str(xi_en)+'.raw.mat'
# sio.savemat( save_file_name )
#
file_name = './data/data.month='+str(id)+'.mat'
load_data = sio.loadmat(file_name)
X = load_data['X']
dmax = load_data['dmax']       ; dmax = int(dmax)
y = load_data['y']             ; y = y.flatten()
y_en = load_data['y_en']       ; y_en = int(y_en)
#
# SET THE PARAMETERS ==========================================================
#
Nsweep = 10**4                         # Number of the mcmc sweeps
Nburnin = 10**3                        # size of the mcmc burn in area
#
# PRIORS ======================================================================
#
a_rho = 1.0
b_rho = 1.0
a_lam = 10.0**(-3)
b_lam = 10.0**(-3)
a_sig2 = 10.0**(-3)
b_sig2 = 10.0**(-3)
#
print( '------------------------') 
print(' PRIORS ')
print( '------------------------') 
print( 'a_rho      : ' , a_rho   ) 
print( 'b_rho      : ' , b_rho   ) 
print( 'a_lam      : ' , a_lam   ) 
print( 'b_lam      : ' , b_lam   ) 
print( 'a_sig2     : ' , a_sig2  ) 
print( 'b_sig2     : ' , b_sig2  ) 
print( '------------------------\n') 
#
# MCMC SEEDS ==================================================================
#
ga_seed = array( dmax*[1] )                            # ge_seed
#
ind = nonzero(ga_seed==1)[0]                                      # beta_seed
beta_seed = array( dmax*[0.0] )                                    # beta_seed
beta_seed[ind] = dot( \
					dot( \
						linalg.inv(dot(X[:,ind].T,X[:,ind])), \
						X[:,ind].T \
						), \
						y.T
					)                                             # beta_seed
#
rho_seed = 0.5                                                      # rho_seed
lam_seed = 0.1                                                      # lambda_seed
sig2_seed = mean((dot(X,beta_seed.T)-y)**2)                    # sig2_seed
#
print( '------------------------') 
print(' MCMC SEEDS ')
print( '------------------------') 
print( 'ga_seed.sum      : ' , sum(ga_seed), "/", ga_seed.shape[0] ) 
print( 'rho_seed         : ' , rho_seed ) 
print( 'sig2_seed        : ' , sig2_seed ) 
print( 'lam_seed         : ' , lam_seed ) 
print( '------------------------') 
error_y = mean((dot(X,beta_seed)-y)**2)
print ('Error of the seed', error_y)
print( '------------------------\n')
#
# MCMC Sampler ================================================================
#
print('MCMC Sampler ... \n')
#
# ... Import module
#
print('... import mcmcmodule  ... \n')
#
import mcmcmodule_L2 as mcmcmodule
#
# ... Initialize the random number generator
#
print('... Initialize the random number generator  ... \n')
#
mcmcmodule.init_genrand( random.randint(99999999) )
for i in range(10) : 
	it = mcmcmodule.rnguniform()
	print(it,) 

print('\n')
#
# ... Declare the arrays
#
print('... Declare the arrays  ... \n')
#
pr_ga = array(dmax*[0.0])            # pr_ga_chain
#
pr_ga_chain = array(Nsweep*[dmax*[0.0]])            # pr_ga_chain
ga_chain = array(Nsweep*[dmax*[0]])            # ga_chain
rho_chain = array(Nsweep*[0.0])            # rho_chain
beta_chain = array(Nsweep*[dmax*[0.0]])        # beta_chain
sig2_chain = array(Nsweep*[0.0])           # sig2_chain
lam_chain = array(Nsweep*[0.0])            # lam_chain
#
# ... Auxiliary values (data) 
#
print('... Auxiliary values (data)  ... \n')
#
# YtY = dot(y,y.T)
# XtY = array( dot(X.T,y.T) )
# XtX = dot(X.T,X)
#
(YtY,XtY,XtX) = mcmcmodule.generate_sufficient_data(X, y, y_en, dmax)
#
YtY = float(YtY)
XtY = XtY.flatten()
#
print( '------------------------') 
print(' AUXILIARY VALUES        ')
print( '------------------------') 
print( 'y_en            : ' , y_en ) 
print( 'YtY             : ' , YtY ) 
print( 'XtY             : ' , XtY ) 
print( 'XtY.rank        : ' , rank(XtY) ) 
print( 'XtY.shape       : ' , shape(XtY) ) 
print( 'XtX             : ' , XtX ) 
print( 'XtX.rank        : ' , rank(XtX) ) 
print( 'XtX.shape       : ' , shape(XtX) ) 
print( '------------------------\n') 
#
# ... Load the seeds
#
print('... Load the seeds  ... \n')
#
ga = ga_seed.copy()
beta = beta_seed.copy()
sig2 = sig2_seed
lam = lam_seed
rho = rho_seed
#
sig2_thr = error_y
#
print( '------------------------') 
print(' BURN-IN SEEDS ')
print( '------------------------') 
print( 'rho         : ' , rho_seed ) 
print( 'sig2        : ' , sig2_seed ) 
print( 'lam         : ' , lam_seed ) 
print( '------------------------\n') 
#
# ... ... Burn-in
#
print ( '... Burn-in ...' )
#
for it in range(Nburnin) :
	#
	if mod(it-1,Nburnin/10.0) == 0 : print (Nburnin-it-1)/(Nburnin/10.0)
	#
	(pr_ga,ga,rho,beta,sig2,lam) \
						= mcmcmodule.mcmcsweep( \
											ga, \
											rho,a_rho,b_rho, \
											beta, \
											sig2,a_sig2,b_sig2, \
											lam,a_lam,b_lam, \
											y_en,YtY,XtY,XtX, \
											0, \
											10000.0, \
											dmax \
	)
	#
	# Record the sample
	#
	pr_ga_chain[it,:dmax] = pr_ga[:dmax]
	ga_chain[it,:dmax] = ga[:dmax]
	rho_chain[it] = rho
	beta_chain[it,:dmax] = beta[:dmax]
	sig2_chain[it] = sig2
	lam_chain[it] = lam				)
#
# ... ... Sample
#
print ( '... Sampling ...' )
#
for it in range(Nsweep) :
	#
	if mod(it-1,Nsweep/10.0) == 0 : print (Nsweep-it-1)/(Nsweep/10.0)
	#
	(pr_ga,ga,rho,beta,sig2,lam) \
					= mcmcmodule.mcmcsweep( \
											ga, \
											rho,a_rho,b_rho, \
											beta, \
											sig2,a_sig2,b_sig2, \
											lam,a_lam,b_lam, \
											y_en,YtY,XtY,XtX, \
											0, \
											10000.0, \
											dmax \
											)
	#
	# Record the sample
	#
	pr_ga_chain[it,:dmax] = pr_ga[:dmax]
	ga_chain[it,:dmax] = ga[:dmax]
	rho_chain[it] = rho
	beta_chain[it,:dmax] = beta[:dmax]
	sig2_chain[it] = sig2
	lam_chain[it] = lam
#
# Save ========================================================================
#
print ( 'Saving...' )
#
# save chain
#
file_name = './results/sample.bma.L2.month='+str(id)+'.mat'
sio.savemat( file_name, \
				{ \
				'n_chain':Nsweep, \
				'dmax':dmax, \
				'pr_ga_chain_bma_L2':pr_ga_chain, \
				'ga_chain_bma_L2':ga_chain, \
				'ga_chain_bma_L2':ga_chain, \
				'beta_chain_bma_L2':beta_chain, \
				'rho_chain_bma_L2':rho_chain, \
				'sig2_chain_bma_L2':sig2_chain, \
				'lam_chain_bma_L2':lam_chain \
				} \
			)
#
# save estimates
#
Pr_ga_bma = mean( pr_ga_chain[1:Nsweep,:], 0 )
Pr_bma = mean( ga_chain[1:Nsweep,:], 0 )
beta_bma = mean( beta_chain[1:Nsweep,:], 0 )
rho_bma = mean( rho_chain[1:Nsweep] )
sig2_bma = mean( sig2_chain[1:Nsweep] )
lam_bma = mean( lam_chain[1:Nsweep] )
#
Pr_ga_bma_median = median( pr_ga_chain[1:Nsweep,:], 0 )
beta_bma_median = median( beta_chain[1:Nsweep,:], 0 )
rho_bma_median = median( rho_chain[1:Nsweep] )
sig2_bma_median = median( sig2_chain[1:Nsweep] )
lam_bma_median = median( lam_chain[1:Nsweep] )
#
file_name = './results/estimates.bma.L2.month='+str(id)+'.mat'
sio.savemat( file_name, \
				{ \
				'Pr_aux_est_bma_L2':Pr_ga_bma, \
				'Pr_est_bma_L2':Pr_bma, \
				'beta_est_bma_L2':beta_bma, \
				'rho_est_bma_L2':rho_bma, \
				'sig2_est_bma_L2':sig2_bma, \
				'lam_est_bma_L2':lam_bma, \
				'Pr_aux_est_bma_median_L2':Pr_ga_bma_median, \
				'beta_est_bma_median_L2':beta_bma_median, \
				'rho_est_bma_median_L2':rho_bma_median, \
				'sig2_est_bma_median_L2':sig2_bma_median, \
				'lam_est_bma_median_L2':lam_bma_median, \
				'beta_seed_bma_L2':beta_seed, \
				'sig2_seed_bma_L2':sig2_seed \
				} \
			)
#
print(' ***END*** ')








