#---------------------------------------------------------------------------
# display the results
#---------------------------------------------------------------------------

# the image

base="Results/$1/fd_jmi"

cat $base"_inverted_image.su" $base"_final_inverted_image.su" | \
suwind key=duse max=1 dt=1 | \
suximage wbox=700 hbox=400 perc=99 legend=1 title="FD data: JMI image - first and last iteration" &

# the velocity model

cat $base"_inverted_velocity.su" $base"_final_inverted_velocity.su" Data/velsm.su | \
suwind key=duse max=1 dt=1 | \
suximage wbox=1000 hbox=400 legend=1 cmap=hsv2 title="FD data: JMI velocity - first and last iteration + True model" &

# misfit function

cat $base"_misfit.su" | \
suxgraph style=normal width=800 height=400 x2beg=0 title="$base misfit function" &

# show the input and modeled data and residual

data=Data/fd_data5.su
res=$base"_residual_imaging.su"

# show residual for FD data input

# find out size of 1 shot

nx=`$DELPHIROOT/bin/convert < $data key=fldr last=1 | surange | head -1 | awk '{print $1}'`

# find out last iteration

duse=`suwind < $res key=fldr max=1 | sugethw duse output=geom | tail -1`

# make selection of shots

fldr="5,12,20,28,35"

# select input

sushw < $data key=fldr a=1 c=1 j=$nx | \
suwind key=fldr max=0 accept=$fldr > tmpdat.su

suwind < $res key=duse max=0 accept=$duse | \
suwind key=fldr max=0 accept=$fldr | \
geom_match geometry=tracf file_ref=tmpdat.su > tmpres.su
