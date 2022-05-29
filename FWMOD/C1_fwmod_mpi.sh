#!/bin/sh

mpirun -np 24 ./fwmod_2d_mpi \
    outfolder=testing \
    velocity=../Results/fugro_model/.su \
    density=Data/den.su \
    source=Data/fwmod_source0.su \
    outfile_label=fwmod \
    multiple_order=5 \
    fmin=100 \
    fmax_upper=600\
    operator_vmin=500 \
    operator_vmax=8000 \
    operator_dv=1 \
    dt=0.0004 \
    dx=2 \
    dz=0.4 \
    size_x=151 \
    size_xtap=50 \
    size_z=201 \
    size_t=256 \
    size_src=151 \
    angle=80 >& log_fwmod

# display some shots

suwind < Data/fwmod_data0.su key=fldr j=30 | \
suximage perc=99.5 wbox=900 title="Some selected shots from FWMod" &
