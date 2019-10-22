#!/bin/bash -l


f2py -m mcmcmodule_L2 -c \
		./mcmcmodule_L2_source/*.f  \
		./mcmcmodule_L2_source/*.f90  \
		--fcompiler=gfortran

for (( MONTH=1; MONTH<=12; MONTH+=1 ))
do
	echo "$MONTH"
	python ./Generate_measurements.py $MONTH > log_Generate_measurements_$MONTH
	python ./Run_MCMC_BMA_L2.py $MONTH > log_Run_MCMC_BMA_L2_$MONTH
	python ./Run_MCMC_MPM_L2.py $MONTH > log_Run_MCMC_MPM_L2_$MONTH
done

python ./Run_Inverse_bma.py > log_Run_Inverse_bma
python ./Run_Inverse_mpm.py > log_Run_Inverse_mpm

python ./Run_Validation_bma.py > Run_Validation_bma
python ./Run_Validation_mpm.py > Run_Validation_mpm





