# %%
import numpy as np, pickle, sys, matplotlib as mpl, pickle, os
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm

FUKA_path = os.getenv('HOME')+'/lib/fuka_dev'
pyFUKA_path = FUKA_path+'/codes/PythonTools'
pyFUKA_libspath = pyFUKA_path+'/lib/'
print(pyFUKA_libspath)
sys.path.append(pyFUKA_libspath)
os.environ["HOME_KADATH)"] = FUKA_path
# %env HOME_KADATH=$FUKA_path

# %%
from fuka_plot_tools.setup_utils import get_reader

# %%
data_base_path="/home/user/lib/falcon_data"
# initial_data_path=pyFUKA_path+'/Example_id/converged_NS_TOTAL_BC.togashi.2.23.-0.4.0.11.dat'
# initial_data_path=FUKA_path+'/codes/reader_test/NS_ISO_UNIFORM_ROT.APR4.madm.1.35.omega.0.0.11.dat'
initial_data_path=data_base_path+"/NS_ISO_UNIFORM_ROT.gam2.madm.1.4.omega.0.0.09.dat"
initial_data_path=data_base_path+"/NS_ISO_UNIFORM_ROT.gam2.madm.1.4.chi.0.2.0.09.dat"
if os.path.isfile(initial_data_path) == True:
    readerISO = get_reader(initial_data_path, ns_iso_uniformrot=True)
else:
    print("{} not found",initial_data_path)

# %%
initial_data_path=data_base_path+"/NS_TOTAL_BC.gam2.madm.1.4.chi.0.0.09.dat"
initial_data_path=data_base_path+"/NS_TOTAL_BC.gam2.madm.1.4.chi.0.2.0.09.dat"
print(initial_data_path)
if os.path.isfile(initial_data_path) == True:
    readerNS = get_reader(initial_data_path, ns=True)
else:
    print("{} not found",initial_data_path)

# %%
readerNS.config

# %%
y1, y2 = [-10, 10]
npts = 256

y_coords = np.linspace(y1, y2, num=npts)

x1, x2 = [-10,10]
x_coords = np.linspace(x1, x2, num=npts)

coords_lst = [[x,y,0] for y in y_coords for x in x_coords]

# %%
varname = "vely"
dataISO = readerISO.getExporterFieldValues(varname, coords_lst)
dataISO = np.array(dataISO)
dataISO = dataISO.reshape(npts,npts)

# %%
dataNS = readerNS.getExporterFieldValues(varname, coords_lst)
dataNS = np.array(dataNS)
dataNS = dataNS.reshape(npts,npts)

# %%
X, Y = np.meshgrid(x_coords, y_coords)

# %%
norm = None #colors.LogNorm(vmin=1e-10, vmax = 1e-2, clip=False) 
#norm = colors.Normalize(vmin=-0.1, vmax = 0.1, clip=False) 

# %%
#print(dataNS)

# %%
plt.rcParams.update({
'figure.figsize'    : [8.0 * 2, 4],
'image.cmap' : 'inferno'
})

# %%
#plt.pcolor(X, Y, dataISO, norm=norm)

# %%


# %%
plt.rcParams.update({
'figure.figsize'    : [8.0 * 2, 4],
'image.cmap' : 'inferno'
})
fig,axs = plt.subplots(1,2)
axs[0].pcolor(X, Y, dataISO, norm=norm)
cs=axs[1].pcolor(X, Y, dataNS, norm=norm)
cbarax,kw = mpl.colorbar.make_axes(axs, location='right', orientation='vertical', aspect=30)
cbar = fig.colorbar(
        cs,
        cax=cbarax, ticklocation='right',
        orientation='vertical', extend='both',
#        ticks=np.linspace(norm.vmin, norm.vmax, num=5)
      )

# %%
plt.show()
quit()
