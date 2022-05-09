import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('rdir', type=str)
parser.add_argument('min_freq', type=float)
args = parser.parse_args()

data1=np.load('Results/%s/jmi_vel-final_vel_model.npy'%(args.rdir))
data=data1[75:226,75:276]*1000
data=np.reshape(data,(151,201))
data.tofile('JMI+FWI-vel.bin')

data2=np.load('Results/%s/vel0-final_vel_model.npy'%(args.rdir))
data=data2[75:226,75:276]*1000
data=np.reshape(data,(151,201))
data.tofile('FWIonly-vel.bin')
