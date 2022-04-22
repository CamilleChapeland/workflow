#!/bin/bash

cd main_2019 
cp Makefile_fwmod Makefile
make remake
cd ../

#pick variables
dx=0.4
dz=0.4
dxrcv=2
dxsrc=2

fmin=100
fmax_lower=120
fmax=600

#pick directory name
dir=fugro_model

#run FWI
echo  ----------------- Running full waveform inversion -----------------
python3 FWI/FWI.py $dir $fmin $fmax $dx $dz

#convert files
echo  ----------------- Converting files to su -----------------
#supaste < 

#run FWM
echo  ----------------- Running full wavefield migration -----------------
./FWM/FD-FWM.sh $mod_dir $dxrcv $dxsrc $fmin $fmax $dz


echo  ----------------- Not running inversion and migration -----------------
exit

