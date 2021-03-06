 -N fdelmod
#PBS -q verylong
#PBS -l nodes=1
#PBS -k eo
#PBS -j eo
#
# Illustrates the effect of (absorbing) tapering of the edges, used in Figure 4 of the manual.

#OK:
dt=0.0002
dx=1
cp=1500
clip=1

makewave file_out=wavelet.su dt=$dt nt=1024 fp=85 shift=1 w=g2 verbose=1

makemod file_base=model.su \
       cp0=$cp ro0=1000 cs0=600 sizex=1000 sizez=1000 \
       dx=$dx dz=$dx orig=0,0 \

export filecp=model_cp.su
export filecs=model_cs.su
export filero=model_ro.su

ntaper=0
../fdelmodc \
        file_cp=$filecp file_cs=$filecs file_den=$filero \
        ischeme=1 \
        file_src=wavelet.su verbose=4 \
        file_rcv=rec.su \
        file_snap=snap_nodisp.su \
        xrcv1=0 xrcv2=1000 dxrcv=5 \
        zrcv1=300 zrcv2=300 \
        rec_type_vx=1 rec_type_pp=1 rec_type_ss=1 rec_int_vx=1 \
        dtrcv=0.004 \
        xsrc=500 zsrc=500 nshot=1 \
        src_type=1 tmod=1.0  \
		ntaper=$ntaper \
		left=3 right=3 bottom=3 top=3 \
        tsnap1=0.2 tsnap2=1.0 dtsnap=0.1 \
        sna_type_ss=1 sna_type_pp=1 fmax=25

suwind key=fldr min=3 max=3 < snap_nodisp_sp.su | \
        supsimage \
        wbox=4 hbox=4 titlesize=-1 labelsize=10 verbose=1 \
        d2=$dx f2=0 clip=$clip \
        label1="depth [m]" label2="lateral position [m]" > snap_tap${ntaper}.eps

supsimage < rec_rp.su \
        wbox=3 hbox=4 titlesize=-1 labelsize=10 clip=$clip verbose=1 \
        label1="time [s]" label2="lateral position [m]" > rec_tap${ntaper}_rp.eps


for ntaper in 50 100;
do
../fdelmodc \
        file_cp=$filecp file_cs=$filecs file_den=$filero \
        ischeme=1 \
        file_src=wavelet.su verbose=4 \
        file_rcv=rec.su \
        file_snap=snap_nodisp.su \
        xrcv1=0 xrcv2=1000 dxrcv=5 \
        zrcv1=300 zrcv2=300 \
        rec_type_vx=1 rec_type_pp=1 rec_type_ss=1 rec_int_vx=1 \
        dtrcv=0.004 \
        xsrc=500 zsrc=500 nshot=1 \
        src_type=1 tmod=1.0  \
		ntaper=$ntaper \
		left=4 right=4 bottom=4 top=4 \
        tsnap1=0.2 tsnap2=1.0 dtsnap=0.1 \
        sna_type_ss=1 sna_type_pp=1 fmax=25

suwind key=fldr min=3 max=3 < snap_nodisp_sp.su | \
        supsimage \
        wbox=4 hbox=4 titlesize=-1 labelsize=10 verbose=1 \
        d2=$dx f2=0 clip=$clip \
        label1="depth [m]" label2="lateral position [m]" > snap_tap${ntaper}.eps

supsimage < rec_rp.su \
        wbox=3 hbox=4 titlesize=-1 labelsize=10 clip=$clip verbose=1 \
        label1="time [s]" label2="lateral position [m]" > rec_tap${ntaper}_rp.eps

done


for npml in 5 10 20;
do
../fdelmodc \
        file_cp=$filecp file_cs=$filecs file_den=$filero \
        ischeme=1 \
        file_src=wavelet.su verbose=4 \
        file_rcv=rec.su \
        file_snap=snap_nodisp.su \
        xrcv1=0 xrcv2=1000 dxrcv=5 \
        zrcv1=300 zrcv2=300 \
        rec_type_vx=1 rec_type_pp=1 rec_type_ss=1 rec_int_vx=1 \
        dtrcv=0.004 \
        xsrc=500 zsrc=500 nshot=1 \
        src_type=1 tmod=1.0  \
		npml=$npml \
		left=2 right=2 bottom=2 top=2 \
        tsnap1=0.2 tsnap2=1.0 dtsnap=0.1 \
        sna_type_ss=1 sna_type_pp=1 fmax=25

suwind key=fldr min=3 max=3 < snap_nodisp_sp.su | \
        supsimage \
        wbox=4 hbox=4 titlesize=-1 labelsize=10 verbose=1 \
        d2=$dx f2=0 clip=$clip \
        label1="depth [m]" label2="lateral position [m]" > snap_pml${npml}.eps

supsimage < rec_rp.su \
        wbox=3 hbox=4 titlesize=-1 labelsize=10 clip=$clip verbose=1 \
        label1="time [s]" label2="lateral position [m]" > rec_pml${npml}_rp.eps

done

cat ../Data/$4/den.su ../Data/$4/truvel.su ../Data/$4/vel0.su | 
exit;

