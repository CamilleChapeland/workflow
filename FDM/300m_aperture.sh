#!/bin/bash 

# define the number of samples, fmin, fmax and dt
ntwav=2750
dtwav=0.00008
fmin=`expr $4 - 50`
fmax=`expr $5 + 50`
t0=0.013
flef=`expr $fmin + 270`
frig=`expr $fmax - 270`

# for modeling
nt=2750
dt=0.00008
xsrc1=0
xsrc2=300
dxsrc=$3
dxrcv=$2
rcv1=0
rcv2=300
minoff=-300
maxoff=300

makewave w=fw fmin=$fmin flef=$flef frig=$frig fmax=$fmax nt=$ntwav dt=$dtwav t0=$t0 > ../$1/wavelet_fdacmod.su

# display the wavelet
suxgraph < ../$1/wavelet_fdacmod.su style=normal title="Wavelet for FD-mod, with shift $t0 seconds" &

cat ../$1/wavelet_fdacmod.su | suspecfx | \
suxgraph style=normal title="Spectrum of wavelet for FD-mod, with shift" &

scale < ../$1/truvel.su a=0 b=1500 | $DELPHIROOT/bin/convert > ../$1/modelhom_cp.su
scale < ../$1/den.su a=0 b=1000 | $DELPHIROOT/bin/convert > ../$1/modelhom_ro.su
# do the FD modeling for a range of shot records
# finally, model the full shot

j=1
/bin/rm ../$1/data.su ../$1/sources.su
xsrc=0

fdelmodc \
   ischeme=1 \
   file_cp=../$1/truvel.su \
   file_den=../$1/den.su \
   file_src=../$1/wavelet_fdacmod.su \
   file_rcv=../$1/shotsfd.su \
   tmod=0.22 \
   xsrc=$xsrc \
   nshot=151 \
   dxshot=$dxsrc \
   zrcv1=0 \
   zrcv2=0 \
   xrcv1=0 \
   xrcv2=300 \
   dxrcv=2 \
   dtrcv=0.00008 \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

fdelmodc \
   ischeme=1 \
   file_cp=../$1/modelhom_cp.su \
   file_den=../$1/modelhom_ro.su \
   file_src=../$1/wavelet_fdacmod.su \
   file_rcv=../$1/shotfdhom.su \
   tmod=0.22 \
   nshot=151 \
   xsrc=$xsrc zsrc=0 \
   zrcv1=0 \
   zrcv2=0 \
   dxshot=$dxsrc \
   xrcv1=0 \
   xrcv2=300 \
   dxrcv=2 \
   dtrcv=0.00008 \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

fdelmodc \
   ischeme=1 \
   file_cp=../$1/modelhom_cp.su \
   file_den=../$1/modelhom_ro.su \
   file_src=../$1/wavelet_fdacmod.su \
   file_rcv=../$1/shotfdhomdir.su \
   tmod=0.22 \
   nshot=151 \
   dxshot=$dxsrc \
   xsrc=$xsrc \
   xrcv1=0 \
   xrcv2=300 \
   dxrcv=2 \
   dtrcv=0.00008 \
   zrcv1=40 \
   zrcv2=40 \
   ntaper=75 tapfact=0.5 \
   taptop=1 tapleft=1 tapright=1 tapbottom=1 \
   verbose=1 

suop2 ../$1/shotsfd_rp.su ../$1/shotfdhom_rp.su| \
rotate trot=-$t0 updatehdr=0| \
pad ntout=$nt  > ../$1/data.su


suwind < ../$1/data.su key=fldr | suximage\

sumute < ../$1/shotfdhomdir_rp.su xmute=-280,0,280 tmute=0.22,0.07,0.22 \
   ntaper=10 mode=1 | sumute xmute=-300,0,300 tmute=0.22,0.015,0.22 \
   ntaper=10  > ../$1/shotfdhomdir.mute.su

taper < ../$1/shotfdhomdir.mute.su nxtaper=10 ntaper=20 | \
pad ntout=$nt | \
taper ntaper=10 nxtaper=10 | \
kxextrap dz=-40 c=1500 |\
rotate trot=-$t0 updatehdr=0 >> ../$1/sources.su

suwind < ../$1/sources.su key=fldr | suximage\

