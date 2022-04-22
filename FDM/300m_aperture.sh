#/bin/bash

# define the number of samples, fmin, fmax and dt
ntwav=3072
dtwav=0.00008
fmin=100
fmax=700
t0=0.026
flef=`expr $fmin + 300`
frig=`expr $fmax - 300`

# for modeling
nt=1250
dt=0.00008
xsrc1=400
xsrc2=700
dxsrc=$3
dxrcv=$2
rcv1=400
rcv2=700
minoff=-300
maxoff=300

inc=`ddiv $dxsrc $dxrcv`
inc=`dnint $inc`

#=========================================================================
# create the time wavelet; shift is needed to make it causal
#=========================================================================
makewave w=fw fmin=$fmin flef=$flef frig=$frig fmax=$fmax nt=$ntwav dt=$dtwav t0=$t0 > $1/wavelet_fdacmod.su

# display the wavelet
suxgraph < $1/wavelet_fdacmod.su style=normal title="Wavelet for FD-mod, with shift $t0 seconds" &

cat $1/wavelet_fdacmod.su | suspecfx | \
suxgraph style=normal title="Spectrum of wavelet for FD-mod, with shift" &

scale < $1/vel.su a=0 b=1500 | $DELPHIROOT/bin/convert > $1/modelhom_cp.su
scale < $1/den.su a=0 b=1000 | $DELPHIROOT/bin/convert > $1/modelhom_ro.su

# do the FD modeling for a range of shot records
# finally, model the full shot
j=1
/bin/rm $1/data.su $1/sources.su

for ((xsrc=$xsrc1;xsrc<$xsrc2;xsrc+=$dxsrc))
do

fdacmod \
   file_vel=$1/vel.su \
   file_den=$1/den.su \
   file_src=$1/wavelet_fdacmod.su \
   file_rcv=$1/shotsfd.su \
   tmod=0.22 \
   xsrc=$xsrc \
   xrcv1=$rcv1 \
   xrcv2=$rcv2 \
   dxrcv=$dxrcv \
   dtrcv=$dt \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

fdacmod \
   file_vel=$1/modelhom_cp.su \
   file_den=$1/modelhom_ro.su \
   file_src=$1/wavelet_fdacmod.su \
   file_rcv=$1/shotfdhom.su \
   tmod=0.22 \
   xsrc=$xsrc zsrc=0\
   xrcv1=$rcv1 \
   xrcv2=$rcv2 \
   dxrcv=$dxrcv \
   dtrcv=$dt \
   zrcv1=0 \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

fdacmod \
   file_vel=$1/modelhom_cp.su \
   file_den=$1/modelhom_ro.su \
   file_src=$1/wavelet_fdacmod.su \
   file_rcv=$1/shotfdhomdir.su \
   tmod=0.22 \
   xsrc=$xsrc zsrc=0\
   xrcv1=$rcv1 \
   xrcv2=$rcv2 \
   dxrcv=$dxrcv \
   dtrcv=$dt \
   zrcv1=40 \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

sumute < $1/shotfdhom.su xmute=-280,0,280 tmute=0.22,0.07,0.22 \
   ntaper=10 mode=1 | sumute xmute=-300,0,300 tmute=0.22,0.015,0.22 \
   ntaper=10 > $1/shotfdhom.mute.su


#suximage < $1/shotfdhom.su perc=99 
#suximage < $1/shotfdhom.mute.su perc=99
suop2 $1/shotsfd.su $1/shotfdhom.su| \
sushw key=fldr a=$j | \
rotate trot=-$t0 updatehdr=0 | \
pad ntout=$nt >> $1/data.su

#suximage < $1/data.su 

sumute < $1/shotfdhomdir.su xmute=-280,0,280 tmute=0.22,0.07,0.22 \
   ntaper=10 mode=1 | sumute xmute=-300,0,300 tmute=0.22,0.015,0.22 \
   ntaper=10  > $1/shotfdhomdir.mute.su

#suximage < $1/shotfdhomdir.mute.su perc=99

taper < $1/shotfdhomdir.mute.su nxtaper=10 ntaper=20 | \
sushw key=fldr a=$j | \
pad ntout=$nt | \
taper ntaper=10 nxtaper=10 | \
kxextrap dz=-40 c=1500 |\
rotate trot=-$t0 updatehdr=0 >> $1/sources.su
#suximage < $1/sources.su

j=$((j+1))

done
