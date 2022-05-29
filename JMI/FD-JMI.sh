#!/bin/bash
mkdir ../logs/$2/

xmin=$7
xmax=$8
dxrcv=$3

xdif=`echo "$xmax-$xmin" | bc`
nrcv=`echo "($xdif/2)+1" | bc`
#nshots=`echo "$nrcv-2" |bc`
nshots=$nrcv
fmin_upper=`echo "$6+50" |bc`

# run the actual JMI
mpirun -n 12 jmi_2d_mpi \
    initial_velocity=../$1/vel02.su \
    data=../$1/data.su \
    source=../$1/sources.su \
    outfile_label=$2/fd_jmi \
    velocity_update_iter_jump=3 \
    reflectivity_const_weight=0.9 \
    output_pershot_info=0 \
    output_residual_info=1 \
    if_model_mask=0 \
    fmin=$6 \
    fmax_lower=$fmin_upper \
    fmax_upper=$5\
    operator_vmin=1000 \
    operator_vmax=4000 \
    operator_dv=1 \
    dt=0.00008 \
    dx=$3 \
    dz=0.4 \
    size_x=$nrcv \
    size_xtap=75 \
    size_z=201 \
    size_t=2500 \
    size_src=$nshots \
    angle=80 \
    size_iter=35,35,35,25,25,25,25,25,25,25\
    data_mask_type=1 \
    velocity_smooth_z=1 \
    velocity_smooth_x=11 \
    velocity_smooth_N=2 \
    illummatrix_smooth_z=3 \
    illummatrix_smooth_x=21 \
    illummatrix_smooth_N=3 >& ../logs/$2/log_fd_jmi
exit
