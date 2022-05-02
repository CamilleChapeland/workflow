#!/bin/sh

# run the actual JMI
mpirun -n 12 jmi_2d_mpi \
    initial_velocity=../$1/vel0.su \
    data=../$1/data.su \
    source=../$1/sources.su \
    outfile_label=fd_jmi \
    velocity_update_iter_jump=10 \
    reflectivity_const_weight=0.9 \
    output_pershot_info=0 \
    output_residual_info=1 \
    if_model_mask=0 \
    fmin=150 \
    fmax_lower=200 \
    fmax_upper=500\
    operator_vmin=1300 \
    operator_vmax=3000 \
    operator_dv=1 \
    dt=0.00008 \
    dx=2 \
    dz=0.4 \
    size_x=151 \
    size_xtap=75 \
    size_z=201 \
    size_t=1500 \
    size_src=151 \
    angle=80 \
    size_iter=50,50,50,50,35,35,35,25,25,25\
    data_mask_type=1 \
    velocity_smooth_z=1 \
    velocity_smooth_x=10 \
    velocity_smooth_N=1 \
    illummatrix_smooth_z=1 \
    illummatrix_smooth_x=3 \
    illummatrix_smooth_N=1 >& ../logs/log_fd_jmi
exit
