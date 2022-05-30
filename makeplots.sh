#!/bin/bash

#which plots would you like to make?
jmi=0
fwi=1
fwm=0
true_model=0
data=0
misfit=0

#Make directory pathways
dir=$1
mod_dir=../Data/$dir
base=Results/$dir

###display JMI results
if [ $jmi -gt 0 ]; then
cat $base/fd_jmi_final_inverted_velocity.su  $base/jmi_vel-final_vel_model.su  $mod_dir/truvel2.su| \
suwind key=duse max=1 dt=1 | \
suximage n1=201 wbox=1200 hbox=400 legend=1 cmap=hsv2 title="Velocity estimate from JMI (left), FWI using JMI starting model (middle) and true velocity model (right)" &

fi

###display FWI results
if [ $fwi -gt 0 ]; then
#cat  $base/jmi_vel-final_vel_model.su $base/vel0-final_vel_model.su | \
#suwind key=duse max=1 dt=1 | \
#suximage wbox=1200 hbox=400 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only and JMI+FWI" &

#cat  $base/vel0-final_vel_model.su $base/jmi_vel-final_vel_model.su $mod_dir/truvel2.su | \
cat $base/vel0-final_vel_model.su $base/jmi_vel-final_vel_model.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=1200 hbox=400 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI, starting velocity model and true velocity model" &

fi

###display FWM results
if [ $fwm -gt 0 ]; then
#make true reflection image
suprod $mod_dir/truvel2.su $mod_dir/den2.su | suop op=refl > $base/refl.su

cat $base/vel0-fwm_final_inverted_image.su  $base/jmi_vel-fwm_final_inverted_image.su |
suximage bclip=0.49 wclip=-0.49 wbox=1200 hbox=400 legend=1 title="Reflectivity using FWI-only background velocity and JMI+FWI background velocity" &

cat $base/refl.su |
suximage bclip=0.49 wclip=-0.49 wbox=1200 hbox=400 legend=1 title="Reflectivity using FWI-only background velocity and JMI+FWI background velocity" &
fi

###display true model
if [ $true_model -gt 0 ]; then
cat $mod_dir/truvel2.su | \
suximage n1=201 wbox=1500 hbox=400 legend=1 cmap=hsv2 title="Velocity" &

cat $base/refl.su | \
suximage n1=201 wbox=1500 hbox=400 legend=1 title="Reflectivity" &
fi

###display misfit functions 
if [ $misfit -gt 0 ]; then
cat $base/fd_jmi_misfit.su $base/vel0-fwm_misfit.su  $base/jmi_vel-fwm_misfit.su| \
suxgraph style=normal width=1200 height=400 x2beg=0 title="Misfit functions for 1. JMI, 2. FWM with vel0 and 3. FWM with jmi_vel" &

cat $base/fd_jmi_misfit.su $base/vel0-fwm_misfit.su  $base/jmi_vel-fwm_misfit.su| \
suxgraph style=normal width=1200 height=400 x2beg=0 title="Misfit functions for 1. JMI, 2. FWM with vel0 and 3. FWM with jmi_vel" &
fi

###display data residuals functions 
if [ $data -gt 0 ]; then

data=$mod_dir/data.su
res=$base/vel0-fwm_residual_imaging.su 

# find out size of 1 shot
nx=`$DELPHIROOT/bin/convert < $data key=fldr last=1 | surange | head -1 | awk '{print $1}'`
# find out last iteration
duse=`suwind < $res key=fldr max=1 | sugethw duse output=geom | tail -1`

# make selection of shots
fldr="75"

#FWI-FWM
sushw < $data key=fldr a=1 c=1 j=$nx | \
suwind key=fldr max=0 accept=$fldr > tmpdat.su

suwind < $res key=duse max=0 accept=$duse | \
suwind key=fldr max=0 accept=$fldr | \
geom_match geometry=tracf file_ref=tmpdat.su > tmpres1.su
suop2 tmpdat.su tmpres1.su > tmpmod.su
n2=`surange < tmpdat.su | head -1 | awk '{print $1}'`

cat tmpdat.su tmpmod.su tmpres1.su | \
suximage width=900 height=600 title="1=$data 2=mod 3=$res"

#JMI-FWI-FWM
res=$base/jmi_vel-fwm_residual_imaging.su 
sushw < $data key=fldr a=1 c=1 j=$nx | \
suwind key=fldr max=0 accept=$fldr > tmpdat.su

suwind < $res key=duse max=0 accept=$duse | \
suwind key=fldr max=0 accept=$fldr | \
geom_match geometry=tracf file_ref=tmpdat.su > tmpres2.su
suop2 tmpdat.su tmpres2.su > tmpmod.su

n2=`surange < tmpdat.su | head -1 | awk '{print $1}'`
cat tmpdat.su tmpmod.su tmpres2.su | \
suximage width=900 height=600 title="1=$data 2=mod 3=$res"

cat tmpres1.su tmpres2.su | \
suximage perc=99.9 width=900 height=600 title="1=$data 2=mod 3=$res"


fi
rm *tmp*
exit

