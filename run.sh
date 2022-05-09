#!/bin/bash

cd ../main_2019 
cp Makefile_fwmod Makefile
make remake
cd ../workflow

#####################################################################
#################### SOME VARIABLES AND FILES #######################
#####################################################################

#pick variables
dx=0.4
dz=0.4
dxrcv=2
dxsrc=2
fmin=150
fmax_lower=200
fmax=600

#pick subsurface model to run
#the model executables must be stored in the ../models/ directory 
#for dir in fugro_model inclusions_1 inclusions_2 karst_1 karst_2 karst_3 water_cave
for dir in fugro_model
do
echo ==================== MODEL $dir =========================

mkdir Results/$dir
mkdir ../Data/$dir

#####################################################################
#################### SELECT WHICH EXECUTABLE ########################
#####################################################################

#pick whether to make synthetic subsurface model
remodel=0

#pick whether to remake FDM data
remakedata=0

#pick whether to run the full workflow
workflow=1

#pick whether to run JMI component
jmi=1
#pick whether to run FWI component
fwi=1
#pick whether to run FWM component
fwm=1

#####################################################################
######################## MODEL AND DATA #############################
#####################################################################

#create subsurface model
if [ $remodel -gt 0 ]; then
../models/$dir'.sh' $dx $dz $dxrcv $dir
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
########################### WORKFLOW  ###############################
#####################################################################

if [ $workflow -gt 0 ]; then

#begin workflow
cd ../main_2019
cp Makefile_jmi Makefile
make remake
cd ../workflow
mkdir ../Results/$dir

############################## JMI ##################################
#run JMI
if [ $jmi -gt 0 ]; then
echo  ----------------- Running joint migration inversion -----------------
./JMI/FD-JMI.sh $mod_dir $dir $dxrcv $dz $fmin $fmax $fmax_lower  
else 
echo ----------------- Re-using JMI results -----------------
fi

cat Results/$dir/fd_jmi_final_inverted_velocity.su  ../$mod_dir/vel02.su ../$mod_dir/truvel2.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=200 legend=1 cmap=hsv2 title="Velocity estimate from JMI (left), starting model (middle) and true velocity model (right)" &

############################## FWI ##################################
#convert files to npy
echo  ----------------- Converting files to bin -----------------
sustrip < Results/$dir/fd_jmi_final_inverted_velocity.su > jmi_vel.bin
sustrip < ../$mod_dir/truvel2.su head=data.head > truvel.bin
sustrip < ../$mod_dir/vel02.su > vel0.bin

#run FWI
if [ $fwi -gt 0 ]; then
echo  ----------------- Running full waveform inversion -----------------
python3 FWI/FWI_smooth.py vel0 $dir $fmin $fmax $dxrcv $dz
python3 FWI/FWI_smooth.py jmi_vel $dir $fmin $fmax $dxrcv $dz
else 
echo  ----------------- Re-using FWI results -----------------
fi

#convert files back to su
echo  ----------------- Converting files to su -----------------
python3 np2su.py $dir $fmin
supaste < JMI+FWI-vel.bin ns=201 head=data.head > Results/$dir/vel0-final_vel_model.su   
supaste < FWIonly-vel.bin ns=201 head=data.head > Results/$dir/jmi_vel-final_vel_model.su   

#display FWI results
cat  Results/$dir/vel0-final_vel_model.su  Results/$dir/jmi_vel-final_vel_model.su ../$mod_dir/vel02.su ../$mod_dir/truvel2.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=200 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI, starting velocity model and true velocity model" &


############################## FWM ##################################
#run FWM
if [ $fwm -gt 0 ]; then
echo  ----------------- Running full wavefield migration -----------------
./FWM/FD-FWM.sh vel0 $dir $dxrcv $dz $fmin $fmax $fmax_lower 
./FWM/FD-FWM.sh jmi_vel $dir $dxrcv $dz $fmin $fmax $fmax_lower 
else 
echo  ----------------- Re-using FWM results -----------------
fi

# display FWM results
base="Results/"
cat $base"vel0-fwm_final_inverted_image" $base"jmi_vel-fwm_final_inverted_image.su" | \
suwind key=duse max=1 dt=1 | \
suximage wbox=800 hbox=400 perc=99 legend=1 title="FD data: FWM first and last iteration" &

else
echo  ----------------- Not running inversion and migration -----------------
./makeplots.sh $dir
fi
rm *.bin
rm *tmp*.su
rm data.head
done
exit

