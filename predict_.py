##########################
## SUPRESS ALL WARNINGS ##
##########################
import warnings, logging, os
warnings.filterwarnings('ignore',category=FutureWarning)
logging.disable(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
##########################

import numpy as np
import subprocess
from subprocess import DEVNULL

# import plotting library
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import sys, getopt
DB_DIR = "/home/krypton/projects/TrR_for_design"

####################
## load libraries ##
####################
sys.path.append(DB_DIR)
from utils import *
from resnet import *
from to_pdb_3 import *

# HACK to fix compatibility issues with rtx2080
config = tf1.ConfigProto()
config.gpu_options.allow_growth = True

## PARSE COMMAND LINE ##
def usage(err=None):
  print("-------------------------------------------------------------------------------------")
  print("TrRosetta predict")
  print("-------------------------------------------------------------------------------------")
  print("-i         :   input fasta")
  print("-o         :   output prefix")
  print("-------------------------------------------------------------------------------------")
  sys.exit(err)

def main(argv):
   fas=None; pre=None
   for opt, arg in getopt.getopt(argv,"i:o:l:n:h",["in=","out="])[0]:
      print(opt,arg)
      if   opt in ("-i","--in"):    fas       = arg
      elif opt in ("-o","--out"):   pre       = arg

   if fas is None or pre is None:
      usage(f"ERROR: I/O not defined")

   headers, seqs = parse_fasta(fas)
   msa = mk_msa(seqs)
   print(f"found {len(seqs)} number of sequences")

   print("Setting up the model")
   design_model = mk_design_model()
   for n,(seq,nam) in enumerate(zip(msa,headers)):
      nam = nam.split(" ")[0]
      out = design_model.get(seq[None])
      seq = N_to_AA(seq.argmax(-1))[0]

      avg_feat = {"theta":[],"phi":[],"cb":[],"omega":[]}

      #save feat
      feats = []
      for o in out:
        feat = split_feat(o["feat"][0])
        feats.append(feat)
        for f in avg_feat:
          avg_feat[f].append(feat[f])

      for f in avg_feat:
        avg_feat[f] = np.mean(avg_feat[f],0)

      feats.append(avg_feat)
      np.save(f"{pre}_{nam}_{n}.npy",feats)

      #save pdb
      xyz = vals_to_xyz(*bins_to_vals(**avg_feat))
      save_PDB(f"{pre}_{nam}_{n}.pdb", xyz, seq)
      subprocess.run(["/home/krypton/projects/TrR_for_design/scwrl4/Scwrl4",
                      "-i",f"{pre}_{nam}_{n}.pdb","-o",f"{pre}_{nam}_{n}.scwrl4.pdb"],stdout=DEVNULL,stderr=DEVNULL)
      print(f"saving: {pre}_{nam}_{n}.pdb")

      #save image
      # plt.figure(figsize=(5*6,5))
      # for k in range(6):
      #   plt.subplot(1,6,k+1)
      #   if k < 5: plt.title(f"model {k}")
      #   else: plt.title(f"avg_model")
      #   plt.imshow(feats[k]["cb"].argmax(-1))
      # plt.savefig(f"{pre}_{nam}_{n}.png", bbox_inches='tight')
      # plt.close()

class mk_design_model:
  '''
  --------------------------------------------------------------------------------
  mk_design_model
  --------------------------------------------------------------------------------
  '''
  ###############################################################################
  # SETUP model to get LOSS and GRADIENTS
  ###############################################################################
  def __init__(self, eps=1e-8):

    K.clear_session()
    K1.set_session(tf1.Session(config=config))
    # inputs
    I = Input(shape=(None, 21), dtype=tf.float32)
    F = RESNET()(I)

    # define model
    self.model = Model(I,F)

    # save weights
    self.weights = []
    for token in ["xaa","xab","xac","xad","xae"]:
      self.weights.append(load_weights(f"{DB_DIR}/models/model_{token}.npy", mode="TrRosetta"))

  def get(self, seq):
    out = []
    for n in range(5):
      self.model.set_weights(self.weights[n])
      feat = self.model.predict(seq)
      cb = feat[...,38:75].argmax(-1)
      out.append({"feat":feat, "cb":cb})
    return out

##############################################################################

def parse_fasta(filename, a3m=False):
  '''function to parse fasta file'''
  if a3m:
    # for a3m files the lowercase letters are removed
    # as these do not align to the query sequence
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))
  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if len(line) > 0:
      if line[0] == ">":
        header.append(line[1:])
        sequence.append([])
      else:
        if a3m: line = line.translate(rm_lc)
        else: line = line.upper()
        sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  return header, sequence

def mk_msa(seqs):
  '''one hot encode msa'''
  alphabet = list("ARNDCQEGHILKMFPSTWYV-")
  states = len(alphabet)
  alpha = np.array(alphabet, dtype='|S1').view(np.uint8)
  msa = []
  for seq in seqs:
    seq = np.array(list(seq), dtype='|S1').view(np.uint8)
    for n in range(states):
       seq[seq == alpha[n]] = n
    seq[seq > states] = states-1
    msa.append(np.eye(states)[seq])
  return msa
##############################################################################
##############################################################################

if __name__ == "__main__":
   main(sys.argv[1:])
