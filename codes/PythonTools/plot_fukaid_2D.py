#!/usr/bin/python3
import numpy as np, pickle, sys, matplotlib as mpl, pickle, os
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm
from matplotlib.gridspec import GridSpec

pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)
from fuka_plot_tools.setup_argparse import *
from fuka_plot_tools.setup_utils import *

def plot_2d(
  x_coords, 
  y_coords, 
  data,
  axs=None,
  plotz=False, **kwargs):
  
  X, Y = np.meshgrid(x_coords, y_coords)
  
  cs = axs.pcolor(X, Y, data, **kwargs)

  axs.set_xlabel(r"$x [{\rm M}_\odot]$",size=LabelSize)
  if plotz:
    axs.set_ylabel(r"$z [{\rm M}_\odot]$",size=LabelSize)
  else:
    axs.set_ylabel(r"$y~[{\rm M}_\odot]$",size=LabelSize)

  axs.set_xticks(ticks=np.linspace(np.min(x_coords),np.max(x_coords), num=5))
  axs.set_yticks(ticks=np.linspace(np.min(y_coords),np.max(y_coords), num=5))
  return cs

def plot_2d_quiver (x_coords,
  y_coords,
  vecx_data, vecy_data,
  axs=None, **kwargs):
  X = np.array(x_coords)
  Y = np.array(y_coords)
  MagnitudeGrid = np.hypot(vecx_data, vecy_data)
  axs.quiver(X, Y, vecx_data, vecy_data,MagnitudeGrid, units='xy',pivot='tail', alpha=1, scale=args.qscale, minlength=1e-3, **kwargs)

if __name__ == "__main__":
  args = get_args()
  f, ispickle = check_ID_filename(args.filename)
  print("Reading from file: {}".format(f))
  
  # setup norm
  if args.vmin is not None and args.vmax is not None:
    norm = colors.Normalize(vmin=args.vmin, vmax = args.vmax, clip=False) 
  else: 
    norm = None

  pltvars  = args.vars if type(args.vars) is list else [args.vars]
  pltqvars = args.qvars if args.qvars is None else [args.qvars]
  
  if not pltqvars is None and len(pltqvars) > 1:
    raise ValueError("Quiver Vars: Quiver only works with one qvar")
  elif not pltqvars is None:
    pltqvars = get_quiver_vars(args.qvars, args.qcomps)

  # Determine if extent is for 1D or 2D plot
  ext = args.extent
  if len(ext) != 4:
    raise ValueError("Extent: incorrect extent information")
  elif len(args.extent) == 4:
    x1, x2, y1, y2 = args.extent
    y_coords = np.linspace(y1, y2, num=args.npts)
    x_coords = np.linspace(x1, x2, num=args.npts)
    qy_coords = np.linspace(y1, y2, num=args.qpts)
    qx_coords = np.linspace(x1, x2, num=args.qpts)

    plotf = lambda data, **kwargs: plot_2d(
      x_coords, 
      y_coords, 
      data,
      norm=norm,
      **kwargs)

  # Setup matplotlib 
  # If a pickle file is used, we assume the 
  # user has setup the plot command in
  # the same way the data was extracted
  plt.close('all')
  height = 4. if args.landscape else 8.
  plt.rcParams.update({
    'figure.figsize'    : [8.0 * len(pltvars), height],
    'image.cmap' : args.cmap
  })
  nrows = 1
  cbarbottom = args.cbar_bottom
  cbarpad = 20
  plt.rcParams['figure.constrained_layout.use'] = True
  fig = plt.figure()
  norm = colors.Normalize(vmin=args.vmin, vmax = args.vmax, clip=False)
  axs = [] # store axes as created
    
  if not ispickle:
    # Setup ID python reader
    reader = get_reader_args(args, f)
    
    for i, v in enumerate(pltvars):    
      inv, sq, plotz, var = parse_var(v)

      print("Extracting data: {} with {}-points".format(var_name(var), args.npts))
      
      # extract data for each var
      if not args.isotropic:
        data = extract_data(
          reader, 
          var, 
          x_coords, 
          y_coords, 
          plotz=plotz,
          square = sq,
          inverse = inv,
          logscale=args.log)
      else:
        data = extract_data_isotropic(
            reader, 
            var, 
            x_coords, 
            y_coords, 
            plotz=plotz,
            square = sq,
            inverse = inv,
            logscale=args.log)
        
      print("Plotting data: {} with {}-shape".format(var_name(var), data.shape))

      # dump data to pickle file
      if args.pickle:
        p, fn = extract_path_filename(f)
        datadim = "2D"
        pdumpf = fn[0:fn.rfind('.')] \
               + "-{}_{}_{}.pickle".format(v,datadim,args.npts)
        with open(pdumpf,'wb') as handle:
          pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Data dumped to {}".format(pdumpf))
      
      cbarlabel = gen_cbarlabel(var_name(var),inv,sq,args.log)
      #print("Plot index: {}".format(1+i))
      ax = fig.add_subplot(nrows, len(pltvars), 1+i)
      
      if len(pltvars) != 1:
        ax.text(
          0.05,
          0.95,
          cbarlabel,
          verticalalignment='top', 
          horizontalalignment='left',
          transform=ax.transAxes,
          fontsize=15,
          color='w')
        cbarlabel = None
      cs = plotf(
        data, 
        axs=ax, 
        plotz=plotz,
      )
      if i > 0:
        ax.set_ylabel("")
        ax.set_yticklabels([])
      ax.set_aspect('auto')
      axs.append(ax)

  else:
    # FIXME this should ideally work for multiple pickle files
    # or the pickle should be able to store multiple datasets
    if len(pltvars) > 1:
      raise ValueError("Too many variables for a pickle file")
    
    inv, sq, plotz, var = parse_var(pltvars[0])
    
    cbarlabel = gen_cbarlabel(var_name(var),inv,sq,args.log)
    with open(f, 'rb') as handle:
      data = pickle.load(handle)
    # store in a list to work in loops later
    axs = [fig.add_subplot(nrows, len(pltvars), 1)]
    cs = plotf(data,axs=axs[0],plotz=plotz)
  
  # Plot colorbar
  if args.cbar:
    font_color = 'w' if args.clean else 'k'
    if len(pltvars) > 1:
      # generate axes across all plots
      loc = 'bottom' if cbarbottom else 'top'
      cbarax,kw = mpl.colorbar.make_axes(axs, location=loc, orientation='horizontal', aspect=70, shrink=1, fraction=0.15,extend='both')

      # https://github.com/matplotlib/matplotlib/issues/22052
      # for horizontal colorbars, we need to invert the axis
      # however there was a bug that caused issues with the
      # extended bars.  See bug report.
      extend = 'both' if mpl.__version__ > '3.5.1' else None
      cbar = fig.colorbar(
        cs,
        cax=cbarax, ticklocation=loc,
        orientation='horizontal', extend='both',
        ticks=np.linspace(norm.vmin, norm.vmax, num=5)
      )
    else:
      cbarax,kw = mpl.colorbar.make_axes(axs, location='right', orientation='vertical', aspect=30)
      cbar = fig.colorbar(
        cs,
        cax=cbarax, ticklocation='right',
        orientation='vertical', extend='both',
        ticks=np.linspace(norm.vmin, norm.vmax, num=5)
      )
      cbar.set_label(cbarlabel, labelpad=cbarpad, rotation=270, color=font_color, size=CbarLabelSize)

  # Plot Quiver
  if not pltqvars is None:
    if args.quiver_filename is None and ispickle:
      raise ValueError("Quiver File: Either set a quiver file -qf <file> or don't use a pickle data file")
    if not args.quiver_filename is None:
      qf, ispickle = check_ID_filename(args.quiver_filename)
    else:
      qf = f
    if not ispickle:
      reader = get_reader_args(args, qf)
      vec_data = []
      for i, v in enumerate(pltqvars):
        inv, sq, plotz, var = parse_var(v)
        print("Extracting Quiver Data: {} with {}-points".format(var_name(var), args.qpts))
        
        # extract data for each var
        data = extract_data(
          reader, 
          var, 
          qx_coords, 
          qy_coords, 
          plotz=plotz)
        vec_data.append(data)
      print("Plotting Quiver data.")
      vecx = vec_data[0]
      vecy = vec_data[1]
      for ax in axs:
        plot_2d_quiver(qx_coords, qy_coords, vecx, vecy,axs=ax, cmap=args.qcmap)
  
  # Misc axes formatting and output
  output_filename = "fig-2D.png"
  if args.clean:
    for ax in axs:
      ax.get_xaxis().set_visible(False)
      ax.get_yaxis().set_visible(False)
      ax.set_axis_off()
    if args.cbar:
      ticks=cbar.get_ticks()
      cbar.set_label(cbarlabel, labelpad=cbarpad, rotation=270, color=font_color, size=CbarLabelSize)

      # set colorbar tick color
      cbar.ax.xaxis.set_tick_params(color=font_color)
      cbar.set_ticklabels(ticks,color=font_color)
      cbar.update_ticks()

      # set colorbar edgecolor 
      cbar.outline.set_edgecolor(font_color)
    plt.savefig(output_filename, dpi=300, bbox_inches="tight", pad_inches=0, facecolor='black')
    print("Plot saved to {}".format(output_filename))
  else:
    for ax in axs:
      ax.minorticks_on()
      ax.tick_params(left=True, bottom=True, top=True, right=True, which='both', labelsize=TickSize)
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print("Plot saved to {}".format(output_filename))
  if args.pltshow:
    plt.show()
