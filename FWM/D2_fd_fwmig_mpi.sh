#!/bin/bash
dxrcv=$2
interpolate < Data/velsm.su d2out=$dxrcv verbose=1 | \
sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/vel_dxrcv.su

mpirun -n 24 jmi_2d_mpi \
    initial_velocity=$1/vel_dxrcv.su \
    data=$1/data.su \
    source=$1/sources.su \
    outfile_label=$1/fd_jmi_inc2_upd \
    velocity_update_iter_jump=0 \
    output_pershot_info=1 \
    output_residual_info=1 \
    if_model_mask=0 \
    fmin=$4 \
    fmax_lower='expr' \
    fmax_upper=$5 \
    operator_vmin=500 \
    operator_vmax=4000 \
    operator_dv=1 \
    dt=0.00008 \
    dx=$dxrcv \
    dz=$6 \
    size_x=151 \
    size_xtap=50 \
    size_z=201 \
    size_t=2200 \
    size_src=151 \
    angle=80 \
    size_iter=10 \
    data_mask_type=1 \
    velocity_smooth_z=1 \
    velocity_smooth_x=1 \
    velocity_smooth_N=1 \
    illummatrix_smooth_z=1 \
    illummatrix_smooth_x=1 \
    illummatrix_smooth_N=1 >& logs/log_fd_fwmig

# display the results

base="Results/fd_jmi_inc2_upd"
cat $base"_inverted_image.su" $base"_final_inverted_image.su" | \
suwind key=duse max=1 dt=1 | \
suximage wbox=700 hbox=400 perc=99 legend=1 title="FD data: FWM first and last iteration" &

# misfit function

cat $base"_misfit.su" | \
suxgraph style=normal width=800 height=400 x2beg=0 title="FD data: FWM misfit function" &
exit
