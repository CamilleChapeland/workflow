#!/bin/bash
# define sampling distance for the modeling
dx=$1
dz=$2

# option for re-modeling or read SEG-Y files
echo "Recreating the model"

# generate a subsurface model for testing JMI code
makemod file_base=model.su cp0=1500 ro0=1000 \
sizex=1000 sizez=80 dx=$dx dz=$dz orig=0,0 \
intt=def poly=0 cp=1550 ro=1700 x=0,450,460,470,480,490,550,600 z=41,41,39,37,36,35,37,45 \
intt=def poly=0 cp=1600 ro=1900 x=550,1000 z=37,37 \
intt=def poly=0 cp=1600 ro=1900 x=400,510,540,550 z=40,41,39,37 \
intt=def poly=0 cp=1900 ro=2000 x=145,205 z=42,42 \
intt=def poly=0 cp=2800 ro=2500 x=190,205,1000 z=48,42,42 \
intt=def poly=0 cp=3000 ro=2200 x=0,145,170 z=42,42,48 \
intt=def poly=0 cp=1600 ro=1700 x=0,470 z=45,45 \
intt=def poly=0 cp=1600 ro=1700 x=490,1000 z=47,47 \
intt=def poly=0 cp=1600 ro=1700 x=470,550 z=45,47 \
intt=def poly=0 cp=1700 ro=1800 x=0,255,300,400,450,500,550,700 z=55,55,54.5,52.5,51,49,47,47 \
intt=def poly=0 cp=1800 ro=1900 x=0,255,300,400,450,500,550 z=60,60,59.5,57.5,56,54,52 \
intt=def poly=0 cp=1900 ro=2000 x=0,255,300,400,450,500,550 z=65,65,64.5,62.5,61,59,57 \
intt=def poly=0 cp=2000 ro=2100 x=0,255,300,400,450,500,550 z=70,70,69.5,67.5,66,64,62 \
intt=def poly=0 cp=1600 ro=1700 x=550,600,1000 z=52,47,47 \
intt=def poly=1 cp=1700 ro=1800 x=500,550,650,750,800,900,966,1000 z=57,52,54,53,52,49,47,47 \
intt=def poly=0 cp=1800 ro=1900 x=450,500,600,700,800,900,1000 z=62,57,59,58,57,54,51 \
intt=def poly=0 cp=1900 ro=2000 x=400,450,550,650,800,900,1000 z=67,62,64,63,62,59,56 \
intt=def poly=0 cp=2000 ro=2100 x=350,400,500,600,800,900,1000 z=72,67,69,68,67,64,61 \
#verbose=1 

mkdir Data/$4

# rename the velocity/demsity file
/bin/mv model_cp.su Data/$4/velf.su
/bin/mv model_ro.su Data/$4/denf.su

# apply smoothing to get starting model
cat Data/$4/velf.su | \
clip minclip=1500 maxclip=2000 | \
smooth ntsm=21 nxsm=31 niter=20 > Data/$4/vel0f.su

suwind min=1001 max=1751 < Data/$4/vel0f.su | sushw key=trwf a=751 | sushw key=tracf,tracl a=1,1 b=1,1 > Data/$4/vel0.su
suwind min=1001 max=1751 < Data/$4/velf.su | sushw key=trwf a=751 | sushw key=tracf,tracl a=1,1 b=1,1 > Data/$4/vel.su
suwind min=1001 max=1751 < Data/$4/denf.su | sushw key=trwf a=751 | sushw key=tracf,tracl a=1,1 b=1,1 > Data/$4/den.su

# also make proper subsampling of the modeling for later use in FWM/JMI
#dxrcv=$2
#interpolate < Data/$4/vel0_5m.su d2out=$dxrcv verbose=1 | \
#sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/$4/vel0f.su
#interpolate < velsm_5m.su d2out=$dxrcv verbose=1 | \
#sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/velsmf.su
#interpolate < Data/$4/vel_5m.su d2out=$dxrcv verbose=1 | \
#sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/$4/velf.su
#interpolate < Data/$4/den_5m.su d2out=$dxrcv verbose=1 | \
#sushw key=gx a=0 b=$dxrcv | sushw key=offset a=0 b=$dxrcv | sushw key=scalco,trid a=0,30 > Data/$4/denf.su

# display the model
cat Data/$4/den.su Data/$4/vel.su Data/$4/vel0.su | 
suximage wbox=1200 hbox=400 title="Density model, Velocity, Smooth velocity and Initial velocity models" legend=1 cmap=hsv2 &
/bin/rm Data/$4/*5m.su
