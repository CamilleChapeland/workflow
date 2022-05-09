#!/bin/bash
#dxrcv=$2
#interpolate < Data/velsm.su d2out=$dxrcv verbose=1 | \
#sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/vel_dxrcv.su

mpirun -n 24 jmi_2d_mpi \
    initial_velocity=Results/$2/$1-final_vel_model.su \
    data=../Data/$2/data.su \
    source=../Data/$2/sources.su \
    outfile_label=$2/$1-fwm \
    velocity_update_iter_jump=0 \
    output_pershot_info=1 \
    output_residual_info=1 \
    if_model_mask=0 \
    fmin=$5 \
    fmax_lower=$7 \
    fmax_upper=$6 \
    operator_vmin=400 \
    operator_vmax=3000 \
    operator_dv=1 \
    dt=0.00008 \
    dx=$3 \
    dz=$4 \
    size_x=151 \
    size_xtap=75 \
    size_z=201 \
    size_t=2750 \
    size_src=151 \
    angle=80 \
    size_iter=10 \
    data_mask_type=1 \
    velocity_smooth_z=1 \
    velocity_smooth_x=1 \
    velocity_smooth_N=1 \
    illummatrix_smooth_z=1 \
    illummatrix_smooth_x=1 \
    illummatrix_smooth_N=1 >& ../logs/$2/log_fd_fwmig

exit
