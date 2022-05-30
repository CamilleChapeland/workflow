#!/bin/bash

cd ../main_2019 
cp Makefile_fwmod Makefile
make remake
cd ../workflow

#####################################################################
#################### SELECT WHICH EXECUTABLE ########################
#####################################################################

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
#pick whether to run FWM component
fwm=1


#####################################################################
#################### SOME VARIABLES AND FILES #######################
#####################################################################

#pick variables
dx=0.4
dz=0.4
xmin=700
xmax=1000
nshots=151
dxrcv=2
dxsrc=2
fmin=150
fmax_lower=200
fmax=600

#pick subsurface model to run
#the model executables must be stored in the ../models/ directory 
for dir in surface_mult
do
echo ==================== MODEL $dir =========================

mkdir Results/$dir
mkdir ../Data/$dir

#####################################################################
######################## MODEL AND DATA #############################
#####################################################################

#create subsurface model
if [ $remodel -gt 0 ]; then
../models/$dir'.sh' $dx $dz $dxrcv $dir $xmin $xmax
echo ----------------- Full subsurface model is generated  -----------------
else 
echo ----------------- Re-using subsurface model -----------------
fi

#model FDM data
mod_dir=Data/$dir
if [ $remakedata -gt 0 ]; then
./FDM/variable_aperture.sh $mod_dir $dxrcv $dxsrc $fmin $fmax $xmin $xmax
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
./JMI/FD-JMI.sh $mod_dir $dir $dxrcv $fmin $fmax $fmax_lower $xmin $xmax  
else 
echo ----------------- Re-using JMI results -----------------
fi

cat Results/$dir/fd_jmi_final_inverted_velocity.su  ../$mod_dir/vel02.su ../$mod_dir/truvel2.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=200 legend=1 cmap=hsv2 title="Velocity estimate from JMI (left), starting model (middle) and true velocity model (right)" &

############################## FWI ##################################
#convert files to npy
echo  ----------------- Converting files to bin -----------------
dx=2
dxrcv=2
interpolate <  Results/$dir/fd_jmi_final_inverted_velocity.su d2out=$dx verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset,f2 a=0,0 b=$dxrcv,$dxrcv | sustrip > Results/$dir/jmi_vel.bin

interpolate <  ../$mod_dir/truvel2.su d2out=$dx verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset,f2 a=0,0 b=$dxrcv,$dxrcv | \
sustrip head=data.head >  Results/$dir/truvel.bin

interpolate <  ../$mod_dir/vel02.su d2out=$dx verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset,f2 a=0,0 b=$dxrcv,$dxrcv | \
sustrip head=data.head >  Results/$dir/vel0.bin


#sustrip < ../$mod_dir/truvel2.su head=data.head >  Results/$dir/truvel.bin
#sustrip < ../$mod_dir/vel02.su >  Results/$dir/vel0.bin
#sustrip < ../$mod_dir/data.su >  Results/$dir/data.bin

#run FWI
if [ $fwi -gt 0 ]; then
echo  ----------------- Running full waveform inversion -----------------
DEVITO_LANGUAGE=openmp python3 FWI/FWI_smooth.py jmi_vel $dir $fmin $fmax $dxrcv $dz $nshots
DEVITO_LANGUAGE=openmp python3 FWI/FWI_smooth.py vel0 $dir $fmin $fmax $dxrcv $dz $nshots
#mpirun -n 1 python3 FWI/FWI_smooth.py vel0 $dir $fmin $fmax $dx $dz $nshots
#mpirun -n 1 python3 FWI/FWI_smooth.py jmi_vel $dir $fmin $fmax $dx $dz $nshots
#python3 FWI/FWI.py jmi_vel $dir $fmin $fmax $dxrcv $dz $nshots
else 
echo  ----------------- Re-using FWI results -----------------
fi

#convert files back to su
echo  ----------------- Converting files to su -----------------
python3 np2su.py $dir $fmin
supaste < JMI+FWI-vel.bin ns=201 head=data.head | \
suwind key=fldr max=149 > Results/$dir/jmi_vel-final_vel_model.su     

supaste < FWIonly-vel.bin ns=201 head=data.head | \
suwind key=fldr max=149 > Results/$dir/vel0-final_vel_model.su 

#display FWI results
cat  Results/$dir/vel0-final_vel_model.su  Results/$dir/jmi_vel-final_vel_model.su ../$mod_dir/truvel2.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=200 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI, starting velocity model and true velocity model" &


############################## FWM ##################################
#run FWM
if [ $fwm -gt 0 ]; then
echo  ----------------- Running full wavefield migration -----------------
#./FWM/FD-FWM.sh vel0 $dir $dxrcv $dz $fmin $fmax $fmax_lower 
./FWM/FD-FWM.sh vel0 $dir $dxrcv $fmin $fmax $fmax_lower $xmin $xmax  
#./FWM/FD-FWM.sh jmi_vel $dir $dxrcv $dz $fmin $fmax $fmax_lower 
./FWM/FD-FWM.sh jmi_vel $dir $dxrcv $fmin $fmax $fmax_lower $xmin $xmax  
else 
echo  ----------------- Re-using FWM results -----------------
fi

# display FWM results
base="Results/"
cat $base"vel0-fwm_final_inverted_image.su" $base"jmi_vel-fwm_final_inverted_image.su" | \
suwind key=duse max=1 dt=1 | \
suximage wbox=800 hbox=400 perc=99 legend=1 title="FD data: FWM first and last iteration" &

else
echo  ----------------- Not running inversion and migration -----------------
#./makeplots.sh $dir
fi
rm *.bin
rm *tmp*.su
rm data.head
done
exit

