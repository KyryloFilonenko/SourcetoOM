#!/bin/sh
#SBATCH --job-name=sourcetoOM
#SBATCH --time=01:00:00
#SBATCH -n 1
#SBATCH --mem 8192M
#SBATCH --licenses sps
#SBATCH --cpus-per-task=1

if [ -f ${THRONG_DIR}/config/supernemo_profile.bash ]; then
	source ${THRONG_DIR}/config/supernemo_profile.bash
fi
snswmgr_load_setup falaise@5.1.5

root making_root.cpp
