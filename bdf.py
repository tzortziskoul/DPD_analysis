import numpy as np
import os
import hoomd
import gsd
from gsd.fl import GSDFile
from hoomd import data
from hoomd import md
import freud
import sys
import gsd.hoomd
import argparse
import matplotlib.pyplot as plt
bond_dist_avg_hist = []
num_bins = 100
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
np.set_printoptions(threshold=sys.maxsize)

def dist_pbc(p1,p2,Box):
  box_dim = Box[0:3]
  box_inv = 1/box_dim
  tilt_xy = Box[3]
  tilt_xz = Box[4]
  tilt_yz = Box[5]
  d = p1 - p2
  images = np.round(d*box_inv)
  d[:,2] -= images[:,2]*box_dim[2]
  d[:,1] -= images[:,2]*box_dim[2]*tilt_yz
  d[:,0] -= images[:,2]*box_dim[2]*tilt_xz
  images = np.round(d*box_inv)
  d[:,1] -= images[:,1]*box_dim[1]
  d[:,0] -= images[:,1]*box_dim[1]*tilt_xy
  images = np.round(d*box_inv)
  d[:,0] -= images[:,0]*box_dim[0]
  return d

parser = argparse.ArgumentParser(description='Calculate bond distribution function from all particles in gsd-file')

non_opt = parser.add_argument_group('mandatory arguments')

non_opt.add_argument('-i','--input', metavar="<gsd file>",  dest='input_file', type=str, required=True,help="input gsd trajectory file")


args = parser.parse_args()

input_file = args.input_file

in_file = gsd.fl.GSDFile(name=args.input_file, mode='r', application='hoomd', schema='hoomd', schema_version=[1,0])
trajectory  = gsd.hoomd.HOOMDTrajectory(in_file)
N_frames = in_file.nframes
indexes = []
frames = []
distances = []
for i,frame in enumerate(trajectory):

    pos = frame.particles.position
    vel = frame.particles.velocity
    ids = frame.particles.typeid
    Box = frame.configuration.box
    bonds = frame.bonds.group

    bond_vecs = dist_pbc(frame.particles.position[bonds[:,1]],frame.particles.position[bonds[:,0]], \
                      frame.configuration.box[0:6])
    bond_distances = np.linalg.norm(bond_vecs,axis=1)

    bond_dist_hist, bin_edges = np.histogram(\
                    bond_distances, bins=num_bins)
    bond_dist_hist = bond_dist_hist / np.sum(bond_dist_hist) # normalizes density to unity
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    shellVols = (4 * np.pi / 3) * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    bond_dist_avg_hist.append(bond_dist_hist / shellVols)


bond_dist_avg_hist = np.array(bond_dist_avg_hist)
bond_dist_avg_hist = np.average(bond_dist_avg_hist,axis=0)
plt.plot(bin_centers,bond_dist_avg_hist)
plt.title('Bond distribution function - Model B')
#plt.savefig('trimmed_bdf_%s.png'%(input_file))
plt.show()
#np.savetxt(args.output_file, np.c_[Final_deb_array],header="# frame, index, distance")
