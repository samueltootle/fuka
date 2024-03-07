import argparse

# FIXME add excision filling extrapolation
# def add_excision_filling_args(parser):
#   # Excision Extrapolation Arguments
#   parser.add_argument(
#     '--interpolation_order', 
#     '--order', 
#     type=int, 
#     help='interpolation order for extrapolation inside excision region', 
#     default='8'
#   )
#   parser.add_argument(
#     '--interpolation_offset', 
#     '--offset', 
#     type=float, 
#     help='relative radial offset interpolation points', 
#     default='0.'
#   )
#   parser.add_argument(
#     '--relative_spacing_dr', 
#     '--dr', 
#     type=float, 
#     help='relative radial spacing of interpolation points', 
#     default='0.4'
#   )
#   parser.add_argument(
#     '--excision', 
#     '--ex', 
#     help='choose excision filling type (spherical, axial)',
#     default=None
#   )

def add_base_arguments(parser):
  parser.add_argument(
    '--filename', 
    "-f", 
    help='Initial data file to plot from', default=None
  )
  parser.add_argument(
    '--clean', 
    action='store_true', 
    help='toggle clean image without axis'
  )
  parser.add_argument(
    '--vars', 
    nargs='*', 
    help="""
    Which var(s) to plot. 
    - A single 1D plot is made or multiple 2D plots depending on a 1D or 2D extent.
    - if a variable is defined with '-z', XZ will be plotted instead of XY.
    - if a variable is defined with inv, the variable will be plotted as 1/var
    - if a variable is defined with sq, the variable will be plotted as var*var
    ex. conf-zsqinv
    """, 
    default=None)
  parser.add_argument('--bhns', action='store_true', help='Use BHNS ID reader')
  parser.add_argument('--bbh' , action='store_true', help='Use BBH ID reader')
  parser.add_argument('--bns' , action='store_true', help='Use BNS ID reader')
  parser.add_argument('--bh'  , action='store_true', help='Use BH ID reader')  
  parser.add_argument('--ns'  , action='store_true', help='Use NS ID reader')
  parser.add_argument('--ns-diffrot'  , action='store_true', help='Use NS Differential Rotation ID reader')
  parser.add_argument('--ns-iso-norot', action='store_true', help='Use NS Isotropic NOROT ID reader')
  parser.add_argument('--ns-iso-uniformrot', action='store_true', help='Use NS Isotropic UNIFORMROT ID reader')
  parser.add_argument('--ns-iso-diffrot', action='store_true', help='Use NS Isotropic DIFFROT ID reader')
  parser.add_argument('--isotropic'  , action='store_true', help='Isotropic solution')
  parser.add_argument('--pickle', 
    action='store_true', 
    help='Disable creation of pickle file after data extraction',
    default=False)

def add_plot_arguments(parser):
  parser.add_argument('--vmax', type=float, help='norm max', default=None)
  parser.add_argument('--vmin', type=float, help='norm min', default=None)
  parser.add_argument('--log', action='store_true', help='toggle log scaling')
  parser.add_argument('--cbar', action='store_true', help='toggle colorbar')
  parser.add_argument(
    '--extent', 
    nargs='*', 
    help='extent: x1 x2 y1 y2', 
    default=None, 
    type=float
  )
  parser.add_argument(
    '--npts', 
    type=int, 
    help='number of interpolation points in x[y/z]',
    default=128
  )
  parser.add_argument('--pltshow', action='store_true', help='toggle plot show')
  parser.add_argument('--landscape', action='store_true', help='set landscape')
  parser.add_argument(
        "--cmap",
        help="Set matplotlib built-in cmap. default: inferno",
        default="inferno",
        type=str,
    )
  parser.add_argument(
    '--cbar_bottom', 
    action='store_true', 
    help='set if colorbar is at the bottom instead of the top. only useful for multiple 2D plots')

  # FIXME - add l2norm
  # parser.add_argument('--L2', action='store_true', help='only compute L2norm')

def add_quiver_plot_arguments(parser):
  parser.add_argument(
    '--quiver_filename', 
    "-qf", 
    help='Initial data file to plot quiver from', default=None
  )
  parser.add_argument(
     '--qpts','--quiver_points', 
    type=int, 
    help='number of interpolation points in x[y/z] for quiver',
    default=32
  )
  parser.add_argument(
    '--qvars','--quiver_vars', 
    nargs='*', 
    help="""
    Which var(s) to use to plot vector fields. 
    - if a qvariable is defined with '-z', XZ components will be used instead of XY.
    """, 
    default=None)
  parser.add_argument(
    '--qscale',
    type=float,
    help='Scale factor for the quiver plot',
    default=0.4)
  parser.add_argument(
    '--qcomps','--quiver_components',
    help='Vector components to plot (e.g. xy, yz)',
    default = 'xy'
  )
  parser.add_argument(
        "--qcmap",
        help="Set quiver cmap to a matplotlib built-in. default: bone",
        default="bone",
        type=str,
    )

def get_args(print_vars=False):
  parser = argparse.ArgumentParser(
    description="Settings for plotting FUKA initial data"
  )
  add_base_arguments(parser)
  # add_excision_filling_args(parser)
  add_plot_arguments(parser)
  add_quiver_plot_arguments(parser)
  
  args = parser.parse_args()
  if not args.isotropic:
    args.isotropic = args.ns_iso_norot or args.ns_iso_diffrot or args.ns_iso_uniformrot
  
  if args.vars == None and print_vars == False:
    raise ValueError("Var(s) must be supplied with --vars <var1 var2 ... varN>")
  
  if args.filename == None:
    raise ValueError("ID filename must be supplied via --filename <> or -f <>")
  
  if args.extent == None and print_vars == False:
    raise ValueError("""
      Extent must be supplied.
      - 1D extent: --extent <x1 x2>. 
      - 2D extent: --extent <x1 x2 y1 y2>"""
    )
  
  if not args.bhns and \
    not args.bbh and \
    not args.bns and \
    not args.bh and \
    not args.ns and \
    not args.ns_diffrot and \
    not args.ns_iso_norot and \
    not args.ns_iso_uniformrot and \
    not args.ns_iso_diffrot:
    raise ValueError("""
      An ID reader must be specified:
        --bbh, --bns, --bh, --ns, etc.  See --help"""
    )
  return args
