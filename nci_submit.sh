#!/usr/bin/env bash
#PBS -P ge3
#PBS -q expressbw
#PBS -l walltime=02:00:00,mem=128GB,ncpus=28,jobfs=100GB
#PBS -l wd
#PBS -j oe

module unload intel-fc intel-cc
module load intel-fc/16.0.3.210
module load intel-cc/16.0.3.210
module load gdal/2.1.3
module load proj/4.8.0
module load geos/3.5.0
module load R/3.4.0

Rscript gwmodel_sirsam.R
