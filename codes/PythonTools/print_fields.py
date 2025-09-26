#!/usr/bin/env python
import sys, os

pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)
from fuka_plot_tools.setup_argparse import *
from fuka_plot_tools.setup_utils import *

if __name__ == "__main__":
  args = get_args(print_vars=True)
  f, ispickle = check_ID_filename(args.filename)
  print("Reading from file: {}".format(f))
  
  # Setup ID python reader
  reader = get_reader_args(args, f)
  
  for v in reader.vars:
      if 'Kadath' in type(reader.vars[v]).__name__:
          print("{:10}: {}".format(v, reader.vars[v]))
