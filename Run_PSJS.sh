#!/bin/bash
sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/EDN_ECP_PSJS_HKY_One_rate_RV_guess_1_nonclock.sh  
sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/EDN_ECP_PSJS_HKY_One_rate_RV_guess_2_nonclock.sh  
sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/EDN_ECP_PSJS_HKY_One_rate_guess_1_nonclock.sh  
sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/EDN_ECP_PSJS_HKY_One_rate_guess_2_nonclock.sh  
