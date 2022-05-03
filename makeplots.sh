#!/bin/bash
dir=fugro_model
mod_dir=../Data/$dir
base=Results/$dir

#display JMI results
cat $base/fd_jmi_final_inverted_velocity.su  $mod_dir/vel0.su $mod_dir/truvel.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=100 legend=1 cmap=hsv2 title="Velocity estimate from JMI (left), starting model (middle) and true velocity model (right)" &

#MISFIT

#display FWI results
cat $base/vel0_final_vel_model.su  $base/jmi_vel_final_vel_model.su $mod_dir/vel0.su $mod_dir/truvel.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=2000 hbox=100 legend=1 cmap=hsv2 title="Velocity estimate comparaison : FWI only, JMI+FWI, starting velocity model and true velocity model" &

# display FWM results
cat $base/vel0-fwm_final inverted_image $base/jmi_vel-fwm_final_inverted_image.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=700 hbox=400 perc=99 legend=1 title="FD data: FWM first and last iteration" &

#display misfit functions 
cat $base/fd_jmi_misfit.su $base/vel0-fwm_misfit.su  $base/jmi_vel_fwm-misfit.su| \
suxgraph style=normal width=1200 height=400 x2beg=0 title="Misfit functions for 1. JMI, 2. FWM with vel0 and 3. FWM with jmi_vel" &
#e
exit

