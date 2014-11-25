#!/bin/bash
clear
rm -rf MatExp
gfortran -I /usr/local/include -I ~/CTEMsoft/Build/src CTEMMatrixExpo.f90 -o MatExp -L /usr/local/lib -L ~/CTEMsoft/Build/Bin -lfortrancl -lCTEMSoftLib -llapack -framework OpenCL
time ./MatExp
