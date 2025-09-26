#!/usr/bin/env python
import numpy as np, pickle, sys, matplotlib as mpl, pickle, os
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm

pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)
from fuka_plot_tools.setup_argparse import *
from fuka_plot_tools.setup_utils import *

def plot_1d(x_coords, var, data, axs=None, fig=None, cbarlabel=None, **kwargs):
  axs.plot(x_coords, data, **kwargs)
  axs.set_xlabel(r"$x [{\rm M}_\odot]$",size=LabelSize)
  if plotz:
    axs.set_ylabel(r"$z [{\rm M}_\odot]$",size=LabelSize)
  else:
    axs.set_ylabel(r"$y~[{\rm M}_\odot]$",size=LabelSize)

if __name__ == "__main__":
  args = get_args()
  f, ispickle = check_ID_filename(args.filename)
  print("Reading from file: {}".format(f))

  # setup norm
  if args.vmin is not None and args.vmax is not None:
    norm = colors.Normalize(vmin=args.vmin, vmax = args.vmax, clip=False)
  else:
    norm = None

  pltvars = args.vars if type(args.vars) is list else [args.vars]

  ext = args.extent
  if len(ext) < 2:
    raise ValueError("Extent: not enough information")

  x1, x2 = ext[0:2]

  x_coords = np.linspace(x1, x2, num=args.npts)
  y_coords = np.array([0.]) if len(ext) == 2 else np.array(ext[2:])

  # Setup matplotlib
  # If a pickle file is used, we assume the
  # user has setup the plot command in
  # the same way the data was extracted
  plt.close('all')
  if args.landscape:
    plt.rcParams.update({
      'figure.figsize'    : [8.0, 4.0]
    })
  fig, axs = plt.subplots(1,1)
  norm = colors.Normalize(vmin=args.vmin, vmax = args.vmax, clip=False)

  if not ispickle:
    # Setup ID python reader
    reader = get_reader_args(args, f)

    for i, v in enumerate(pltvars):
      inv, sq, plotz, var = parse_var(v)

      print("Plotting: {} with {}-points".format(var_name(var), args.npts))

      for y in y_coords:
        # extract data for each var
        if not args.isotropic:
          data = extract_data(
            reader,
            var,
            x_coords,
            [y],
            plotz=plotz,
            square = sq,
            inverse = inv,
            logscale=args.log,
            absolute=args.abs)
        else:
          data = extract_data_isotropic(
            reader,
            var,
            x_coords,
            [y],
            plotz=plotz,
            square = sq,
            inverse = inv,
            logscale=args.log,
            absolute=args.abs)
        print("Plotting data with shape {}".format(data.shape))

        # dump data to pickle file
        if args.pickle:
          p, fn = extract_path_filename(f)
          datadim = "1D"
          pdumpf = fn[0:fn.rfind('.')] + "-{}_{}.pickle".format(v,datadim)
          with open(pdumpf,'wb') as handle:
            pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
          print("Data dumped to {}".format(pdumpf))
        cbarlabel = gen_cbarlabel(var_name(var),inv,sq,args.log)
        axs.plot(x_coords, data,label=cbarlabel)

  else:
    # FIXME this should ideally work for multiple pickle files
    # or the pickle should be able to store multiple datasets
    if len(pltvars) > 1:
      raise ValueError("Too many variables for a pickle file")

    inv, sq, plotz, var = parse_var(pltvars[0])
    print(var_name(var))
    cbarlabel = gen_cbarlabel(var_name(var),inv,sq,args.log)
    with open(f, 'rb') as handle:
      data = pickle.load(handle)
    cbarlabel = gen_cbarlabel(var_name(var),inv,sq,args.log)
    axs.plot(x_coords, data,label=cbarlabel)

  axs.set_xlabel(r"$x [{\rm M}_\odot]$",size=LabelSize)
  axs.set_ylabel(r"$\mathcal{X}$",size=LabelSize)
  axs.legend(fontsize=LabelSize)

  if args.clean:
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gca().set_axis_off()
    plt.savefig("fig-1D.png", dpi=300, bbox_inches="tight", pad_inches=0, facecolor='black')
  else:
    axs.minorticks_on()
    axs.tick_params(left=True, bottom=True, top=True, right=True, which='both', labelsize=TickSize)
    plt.savefig("fig-1D.png", dpi=300, bbox_inches="tight")
  if args.pltshow:
    plt.show()
