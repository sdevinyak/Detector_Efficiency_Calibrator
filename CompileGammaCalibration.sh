########################################
#Author: Sophia Devinyak (s.devinyak@gmail.com , sdevinya@uwaterloo.ca)
#Date created: November 2021 at TRIUMF under the supervision of Dr. Roger Caballero-Folch



#!/bin/bash

echo Removing old files
rm input.cxx 2>/dev/null
rm input_rdict.pcm 2>/dev/null

echo Compiling the program
rootcint input.cxx -c input.h
g++ `root-config --cflags --glibs` -O3 -lCore -lRIO -lNet -lHist -lGraf -lGui -lSpectrum -DNOROOT=1 -Wunused-variable -Wall -Wextra input.cxx DetectorEfficiency.C -o EfficiencyCalibrator


echo The program is ready to use

##chmod +x CompileGammaCalibration.sh in the terminal to make it excecutable
