import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
import argparse, gsd, gsd.hoomd
import sys
import math

plt.rcParams.update({
    'mathtext.fontset': 'cm',
    'font.family': 'STIXGeneral',
    'xtick.labelsize': 9,
    'text.usetex': True,
    'savefig.bbox': "tight",
    'savefig.pad_inches': 0.1,
    'legend.frameon': False,
    'legend.shadow': False
})

def dist_pbc(x1,x0,Box):
    delta = np.abs(x1 - x0)
    delta= np.where(delta > 0.5 * Box, Box - delta, delta)
    return np.sqrt(np.sum(delta**2.0))

def dist_vec_pbc(x1,x0,Box):
    delta = x1 - x0
    delta= np.where(delta > 0.5 * Box, Box - delta, delta)
    delta= np.where(delta <- 0.5 * Box, Box + delta, delta)
    return delta

def unit_vector(vector):
    if np.linalg.norm(vector)>0:
        return vector / np.linalg.norm(vector)
    else:
        return np.array([0,0,0])

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

def calc_adf(pos,nbins):
    g = []
    for k in np.linspace(0, 180, nbins+1 ) [:-1]:

        dk = 180 / nbins
        total = 0
        for j in range(len(pos)-2):
            v1=dist_vec_pbc(pos[j],pos[j+1],Box)
            v2=dist_vec_pbc(pos[j+2],pos[j+1],Box)
            angle=angle_between(v1,v2)
            all_angle.append(angle)
            total += ((k <= angle) &  (angle < k + dk)).sum()
        g.append(total)
    return np.array(g)


############################### PARSE ##################################

parser = argparse.ArgumentParser(description='Calculate clusters on particles in gsd-file')

non_opt = parser.add_argument_group('mandatory arguments')

non_opt.add_argument('-i','--in',metavar='<gsd-file>', dest='input_file',type=str,\
                required=True,help='input trajectory *.gsd file')
args = parser.parse_args()
input_file = args.input_file
nbins = 180

all_angle =[]
# read gsd file
in_file = gsd.fl.GSDFile(name=input_file, mode='r', application='hoomd', schema='hoomd', schema_version=[1,0])
trajectory  = gsd.hoomd.HOOMDTrajectory(in_file)
N_frames = int(in_file.nframes)
a = []
for i,frame in enumerate(trajectory[-1:]):

       Box = frame.configuration.box[0:3]
       pos = frame.particles.position
       vel = frame.particles.velocity
       ids = frame.particles.typeid
       adf = calc_adf(pos,nbins)
x = np.linspace(0,180,nbins)
#spl = make_interp_spline(x, power, k=3)  # type: BSpline
plt.plot(x,adf)
plt.title('Bond distribution function - Model B')

#plt.hist(adf, bins = np.linspace(0, 180, nbins ))
plt.show()
