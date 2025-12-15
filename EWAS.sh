#!/bin/bash
#SBATCH -A Research_Project-T127716 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p sq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=jp789@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err




beta_file=./Betas/IPSC.Betas_CV_vs_PLCG2_KO.rds                   #rds file contains beta matrix
pheno_file=./Phenotypes/IPSC.Pheno_CV_vs_PLCG2_KO.csv             #phenotype csv file contains samples information 
trait=Genotype                                                    # target column name in phenotype csv file
covars_fact=Plate.Location                                        # Factor covariates to in the regression model
covars_num=""                                                     # Numeric covariates in the regression model                                              
model_lm=~Genotype+Plate.Location                                 # Linear regression model fo running the EWAS analysis
num_threads=16						                              # number of thread to run EWAS analysis in parallel mode 
out_pref=./Jack_EWAS                                              #Out put files prefix

script_dir=.                                                      # R scripts directory 

Rscript ${script_dir}/EWAS.R "$beta_file" "$pheno_file" "$trait" "$covars_fact" "$covars_num" "$model_lm" "$num_threads" "$out_pref" 


