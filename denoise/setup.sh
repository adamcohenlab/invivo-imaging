#!/bin/bash

module load Anaconda3/5.0.1-fasrc02
module load intel/17.0.4-fasrc01
source activate invivo

export TREFIDE="$PWD/../../../Programs/trefide"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TREFIDE/src"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TREFIDE/src/proxtv"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TREFIDE/src/glmgen/lib"

export INCLUDE=${MKL_HOME}/include:$INCLUDE
export LIBRARY_PATH=${MKL_HOME}/lib/intel64:$LIBRARY_PATH
export CPATH=${MKL_HOME}/include:$CPATH
export LD_PRELOAD=${MKL_HOME}/lib/intel64/libmkl_core.so:${MKL_HOME}/lib/intel64/libmkl_sequential.so


