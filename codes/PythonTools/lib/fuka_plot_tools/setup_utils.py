import sys, os
import matplotlib.pyplot as plt, matplotlib
pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)

LabelSize=15
TickSize=15
CbarLabelSize=15

# Set default matplotlib settings
plt.rcParams.update({
    'figure.figsize'    : [8.0, 8.0],
    'text.usetex'       : matplotlib.checkdep_usetex(True),
    'font.family'       : "sans-serif",
    'font.serif'        : "cm",
    'xtick.major.size'  : 6,
    'xtick.minor.size'  : 3,
    'ytick.major.size'  : 6,
    'ytick.minor.size'  : 3,
    'xtick.top' : True,
    'ytick.right' : True,
    'axes.labelsize' : LabelSize,
    'legend.fontsize' : LabelSize * 0.75,
})


def extract_path_filename(f):
  # extract filename and path
  idx = f.rfind('/')
  if idx != -1:
    path = f[0:idx+1]
    f = f[idx+1:]
  else:
    path = './'
  return path, f

def gen_cbarlabel(var_name,sq=False,inv=False,log=False):
  cbarlabel = var_name
  if sq:
    cbarlabel += r'^2'
  if inv:
    cbarlabel = r'1 / ' + cbarlabel
  if log:
    cbarlabel = r'\log \left( '+cbarlabel+r' \right)'
  cbarlabel = r'$'+cbarlabel+r'$'
  return cbarlabel

def check_ID_filename(input_filename):
  path, f = extract_path_filename(input_filename)
  
  # make sure we have the correct ext
  ext = f.split('.')[-1].lower()
  ispickle = (ext == 'pickle')
  
  if ext != 'dat' and not ispickle:
    ext = '.' + f.split('.')[-1]
    f = f.strip(ext) + '.dat'
  
  return path+f, ispickle

def get_reader_args(args, filepathabs):
  return get_reader(filepathabs,
    bhns=args.bhns,
    bbh = args.bbh,
    bns = args.bns,
    ns = args.ns,
    ns_iso_norot = args.ns_iso_norot,
    ns_iso_diffrot = args.ns_iso_diffrot,
    ns_iso_uniformrot = args.ns_iso_uniformrot,
    bh = args.bh,
    )
  
def get_reader(filepathabs,
               bhns=False,
               bns=False,
               bbh=False,
               ns=False,
               ns_iso_norot=False,
               ns_iso_uniformrot=False,
               ns_iso_diffrot=False,
               bh=False):
  reader = None
  if bhns:
    from fukaID_readers.bhns import bhns_reader
    reader = bhns_reader(filepathabs)
  elif bbh:
    from fukaID_readers.bbh import bbh_reader
    reader = bbh_reader(filepathabs)
  elif bh:
    from fukaID_readers.bh import bh_reader
    reader = bh_reader(filepathabs)
  elif bns:
    from fukaID_readers.bns import bns_reader
    reader = bns_reader(filepathabs)
  elif ns:
    from fukaID_readers.ns import ns_reader
    reader = ns_reader(filepathabs)
  elif ns_iso_norot:
    from fukaID_readers.ns import ns_isotropic_norot_reader
    reader = ns_isotropic_norot_reader(filepathabs)
  elif ns_iso_uniformrot:
    from fukaID_readers.ns import ns_isotropic_uniformrot_reader
    reader = ns_isotropic_uniformrot_reader(filepathabs)
  elif ns_iso_diffrot:
    from fukaID_readers.ns import ns_isotropic_diffrot_reader
    reader = ns_isotropic_diffrot_reader(filepathabs)
    
  return reader

def parse_var(var):
  inv = False
  sq =  False
  plotz=False
  if var[-2:] == "sq":
    sq = True
    var = var[:-2]
  if var[-3:] == "inv":
    inv = True
    var = var[:-3]
  if var[-2:] == "-z":
    plotz=True
    var = var[:-2]
  return inv, sq, plotz, var

name_dict = {
  'conf' : r'\Psi',
  'lapse' : r'\alpha',
  'shift' : r'\beta',
  'cPsi' : r'\mathcal{C}_{\Psi}',
  'cLapsePsi' : r'\mathcal{C}_{\alpha\Psi}',
  'cShift_1' : r'\mathcal{C}_{\beta}^1',
  'cShift_2' : r'\mathcal{C}_{\beta}^2',
  'cShift_3' : r'\mathcal{C}_{\beta}^3',
  'logh' : r'\log \left( h \right)',
  'drPsi' : r'\partial_r \Psi',
  'ddrPsi' : r'\partial_r^2 \Psi',
  'A' : r'\hat{A}^i_i',
  'A_11' : r'\hat{A}_{11}',
  'A_12' : r'\hat{A}_{12}',
  'A_13' : r'\hat{A}_{13}',
  'A_22' : r'\hat{A}_{22}',
  'A_23' : r'\hat{A}_{23}',
  'A_33' : r'\hat{A}_{33}',
  'A_xx' : r'\hat{A}_{ij} e^x e^x',
  'A_yy' : r'\hat{A}_{ij} e^y e^y',
  'A_zz' : r'\hat{A}_{ij} e^z e^z',
  'A_xy' : r'\hat{A}_{ij} e^x e^y',
  'A_yz' : r'\hat{A}_{ij} e^y e^z',
  'A_xz' : r'\hat{A}_{ij} e^x e^z',
  'rho' : r'\rho',
  'eps' : r'\epsilon',
  'press' : 'P',
  'dHdx' : 'e^x D_x \log \left( H \right)',
  'P/rho' : r'\frac{P}{\rho}',
  'Stilde' : r'\tilde{S}',
  'Etilde' : r'\tilde{E}',
  'vel_1' : r'U^1',
  'vel_2' : r'U^2',
  'vel_3' : r'U^3',
  'shift_1' : r'\beta_1',
  'shift_2' : r'\beta_2',
  'shift_3' : r'\beta_3',
}

def var_name(var):
  if var in name_dict.keys():
    return name_dict[var]
  else: 
    return var

def extract_data(
  reader, 
  var, 
  x_coords,
  y_coords,
  plotz=False,
  logscale=False,
  square=False,
  inverse=False,
  zval=0):
  
  import numpy as np

  xlen = len(x_coords)
  ylen = len(y_coords)
  
  if not plotz:
    coords_lst = [[x, y, zval] for y in y_coords for x in x_coords]
  else:
    coords_lst = [[x, zval, z] for z in y_coords for x in x_coords]
  
  # FIXME add excision filling extrapolation
  # if excision_extraction is not None:
    # if excision_extraction == "spherical":
    #   data = bh.getExFieldVal_sphere_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)
    # elif excision_extraction == "oldaxial":
    #   data = bh.getExFieldVal_axial_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)
    # elif excision_extraction == "axial":
    #   data = bh.getExFieldVal_spectral_axial_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)

    # else:
  data = reader.getFieldValues(var, coords_lst, -1)
  data = np.array(data)
  
  if square:
    data *= data

  if inverse:
    data = 1. / data

  if logscale:
    data = np.array(np.log10(np.abs(data)+1e-15))
  
  if xlen == ylen:
    data=data.reshape(xlen,ylen)
  return data

def extract_data_isotropic(
  reader, 
  var, 
  x_coords,
  y_coords,
  plotz=False,
  logscale=False,
  square=False,
  inverse=False,
  zval=0):
  
  import numpy as np

  xlen = len(x_coords)
  ylen = len(y_coords)
  
  if not plotz:
    coords_lst = [[x, y] for y in y_coords for x in x_coords]
  else:
    coords_lst = [[x, zval] for z in y_coords for x in x_coords]

  data = reader.getFieldValues(var, coords_lst, -1)
  data = np.array(data)
  
  if square:
    data *= data

  if inverse:
    data = 1. / data

  if logscale:
    data = np.array(np.log10(np.abs(data)+1e-15))
  
  if xlen == ylen:
    data=data.reshape(xlen,ylen)
  return data

def get_quiver_vars(var, plane):
  qvars = list()
  for c in plane:
    if c.lower() == "x" or c.lower() == '1':
      qvars.append(var[0]+"_1")
    elif c.lower() == "y" or c.lower() == '2':
      qvars.append(var[0]+"_2")
    elif c.lower() == "z" or c.lower() == '3':
      qvars.append(var[0]+"_3")
    else:
       ValueError("Quiver Components: Quiver only works with xyz components")
  if len(qvars) > 2:
    ValueError("Quiver Components: Quiver only works with 2D components")
  return qvars
