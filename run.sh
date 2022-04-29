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

fmin=150
fmax_lower=200
fmax=600

#pick directory name
dir=fugro_model

#pick whether to make synthetic subsurface model
remodel=1

#pick whether to remake FDM data
remakedata=1

#pick whether to run the full workflow
workflow=1

#pick whether to run JMI component
jmi=1
#pick whether to run FWI component
fwi=1
#pick whether to run FDM component
fwm=0




#####################################################################
#####################################################################
#make dataset if necessary
#create subsurface model
if [ $remodel -gt 0 ]; then
./models/$dir'.sh' $dx $dz $dxrcv $dir
echo ----------------- Full subsurface model is generated  -----------------
else 
echo ----------------- Re-using subsurface model -----------------
fi

#model FDM data
mod_dir=Data/$dir
if [ $remakedata -gt 0 ]; then
./FDM/300m_aperture.sh $mod_dir $dxrcv $dxsrc $fmin $fmax
echo ----------------- FDM data created -----------------
else 
echo ----------------- Re-using FDM data -----------------
fi

#####################################################################
#####################################################################
if [ $workflow -gt 0 ]; then

#begin workflow
cd main_2019
cp Makefile_jmi Makefile
make remake
cd ../
mkdir Results/$dir

#run JMI
if [ $jmi -gt 0 ]; then
echo  ----------------- Running joint migration inversion -----------------
./JMI/FD-JMI.sh $dir $dxrcv $dxsrc $fmin $fmax $fmax_lower $dz 
./JMI/plot_results.sh $dir
else 
cat  Results/$dir/fd_jmi_final_inverted_velocity.su  $mod_dir/vel0.su $mod_dir/truvel.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=400 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI,tTrue vel, starting vel" &
echo ----------------- Re-using JMI results -----------------
fi

#convert files
echo  ----------------- Converting files to bin -----------------
sustrip < Results/$dir/fd_jmi_final_inverted_velocity.su > jmi_vel.bin
sustrip < $mod_dir/truvel.su head=data.head > truvel.bin
sustrip < $mod_dir/vel0hom.su > vel0.bin

#run FWI
if [ $fwi -gt 0 ]; then
echo  ----------------- Running full waveform inversion -----------------
python3 FWI/FWI.py $dir $fmin $fmax $dxrcv $dz 
python3 FWI/FWI_smooth.py $dir $fmin $fmax $dxrcv $dz
else 
echo  ----------------- Re-using FWI results -----------------
fi
#convert files
echo  ----------------- Converting files to su -----------------
python3 np2su.py $dir $fmin
supaste < JMI+FWI-vel.bin ns=201 head=data.head > Results/$dir/$fmin-JMI+FWI-vel.su   
supaste < FWIonly-vel.bin ns=201 head=data.head > Results/$dir/$fmin-FWIonly-vel.su   

cat  Results/$dir/$fmin-FWIonly-vel.su  Results/$dir/$fmin-JMI+FWI-vel.su $mod_dir/vel0.su $mod_dir/truvel.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=400 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI,tTrue vel, starting vel" &


#run FWM
if [ $fwm -gt 0 ]; then
echo  ----------------- Running full wavefield migration -----------------
./FWM/FD-FWM.sh $mod_dir $dxrcv $dxsrc $fmin $fmax $dz
else 
echo  ----------------- Re-using FWM results -----------------
fi

else
echo  ----------------- Not running inversion and migration -----------------
fi
rm *.bin
exit

