#!/bin/sh

n1out=`suwind < Data/$1/vel0.su count=1 | sugethw ns output=geom`
n2out=`surange  < Data/$1/vel0.su | head -1 | awk '{print $1}'`
dxrcv=$2
dz=$7
dx=$dxrcv

suwind min=1001 max=1751 < Data/$1/vel0f.su | sushw key=trwf a=751 | sushw key=tracf,tracl a=1,1 b=1,1 | \
interpolate d2out=$dxrcv verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/$1/vel0.su

suwind min=1001 max=1751 < Data/$1/velf.su | sushw key=trwf a=751 | sushw key=tracf,tracl a=1,1 b=1,1 | \
interpolate d2out=$dxrcv verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/$1/truvel.su

# run the actual JMI
mpirun -n 12 jmi_2d_mpi \
    initial_velocity=Data/$1/vel0.su \
    data=Data/$1/data.su \
    source=Data/$1/sources.su \
    outfile_label=$1/fd_jmi \
    velocity_update_iter_jump=2 \
    reflectivity_const_weight=0.9 \
    output_pershot_info=0 \
    output_residual_info=1 \
    if_model_mask=0 \
    fmin=$4 \
    fmax_lower=$6 \
    fmax_upper=$5\
    operator_vmin=1300 \
    operator_vmax=3500 \
    operator_dv=1 \
    dt=0.00008 \
    dx=$dx \
    dz=$dz \
    size_x=151 \
    size_xtap=75 \
    size_z=201 \
    size_t=2200 \
    size_src=151 \
    angle=80 \
    size_iter=25,25,25,25,25,25,25,25,25,25\
    data_mask_type=1 \
    velocity_smooth_z=1 \
    velocity_smooth_x=1 \
    velocity_smooth_N=1 \
    illummatrix_smooth_z=3 \
    illummatrix_smooth_x=3 \
    illummatrix_smooth_N=3 >& logs/log_fd_jmi
exit
