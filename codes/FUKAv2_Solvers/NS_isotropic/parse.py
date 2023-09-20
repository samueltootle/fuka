#!/usr/bin/python3
import os, sys, glob
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm
pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)
# from fuka_plot_tools.setup_argparse import *
from fuka_plot_tools.setup_utils import *
import math

import glob

if __name__ == "__main__":
  print("test", glob.glob("res_seq*"))
  reslist=[13] #, 17, 21, 25, 33]
  for R in reslist:
    for d in sorted(glob.glob("seq_rr*A-1")):
      fabs=d+"/NS_ISO_DIFF_ROT.gam2.keh.nc.*.info"
      print(fabs)
      files = glob.glob(fabs)
  #    print(files)
      rho = []
      Mb = []
      res = []
      res_r = []
      Madm = []
      Mk = []
      for f in sorted(files):
        fID, p = check_ID_filename(f)
        ns = get_reader(fID, ns_iso_diffrot=True)
        # print(ns.config['ns']['res'])
        if not ns.config['ns']['res'] == int(R) :
          continue
        # print(ns.vars['nc'], ns.vars['Mb'])
        rho.append(ns.vars['nc']) # * 6.17714e17)
        Mb.append(ns.vars['Mb'])
        Madm.append(ns.vars['Madm'])
        Mk.append(ns.vars['Mk'])
        res_sum = 0
        for reslist in ns.vars['domain_resolutions']:
          res_sum += reslist[0] * reslist[1]
        res.append(res_sum)
        res_r.append(ns.config['ns']['res'])
      plt.plot(rho, Madm, label=d)
      # plt.plot([r**(1./2.) for r in res[:-1]], [math.fabs(1. - M / Madm[-1]) for M in Madm[:-1]], label=d)
      # plt.plot([r**(1./2.) for r in res], [math.fabs(1. - mk / madm) for mk, madm in zip(Mk, Madm)], label=d)
      print(len(rho))
  plt.xlabel(r'$\rho {\rm [g / {cm}^3]}$')
  plt.ylabel(r'$M_b {\rm [M_\odot]}$')
  plt.ylim(-0.1, 2.6)
  plt.legend()
  # plt.xlabel(r'$\bar{N}^{1/2}$')
  # plt.ylabel(r'$|1 - M^{\bar{N}}_{\rm ADM} / M^{\bar{N}, HR}_{\rm ADM}|$')
  # plt.yscale('log')
  plt.minorticks_on()
  plt.savefig('seq1.pdf')