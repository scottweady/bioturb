#!/bin/bash

# this script reads jobsub.slurm and submit it to slurm
# sets the folder name as the slurm job name
foldername=${PWD##*/} # some bash magic
rm ./*.subscript
cp ./jobsub.slurm ./$foldername.subscript
sbatch -p ccb --qos=ccb ./$foldername.subscript
