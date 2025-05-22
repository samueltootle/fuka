#!/usr/bin/python3
import numpy as np, pickle, sys, matplotlib as mpl, pickle, os
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm

pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)
from fuka_plot_tools.setup_utils import get_reader
def test_bns():
  from fukaID_readers.bns import bns_reader

  bns = bns_reader('Example_id/converged_BNS_TOTAL.togashi.30.6.0.0.2.8.q1.0.0.09.dat')
  return bns

def test_bhns():
  from fukaID_readers.bhns import bhns_reader

  bhns = bhns_reader('Example_id/converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11.dat')
  return bhns

def test_bbh():
  from fukaID_readers.bbh import bbh_reader

  bbh = bbh_reader('Example_id/converged_BBH_TOTAL_BC.10.0.0.1.q1.0.0.09.dat')
  return bbh

def test_bh():
  from fukaID_readers.bh import bh_reader

  bh = bh_reader('Example_id/converged_BH_TOTAL_BC.0.5.0.0.09.dat')
  return bh

def test_ns():
  from fukaID_readers.ns import ns_reader

  ns = ns_reader('Example_id/converged_NS_TOTAL_BC.togashi.2.23.-0.4.0.11.dat')
  return ns

f="<ID>.dat"
#readerISO = get_reader(f, ns_iso_diffrot=True, ns_iso_uniformrot=False)
#print(readerISO.vars)
