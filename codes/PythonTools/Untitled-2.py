# %%
import numpy as np, pickle, sys, matplotlib as mpl, pickle, os
import matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib import ticker, cm

FUKA_path = os.getenv('HOME')+'/lib/fuka_dev'
pyFUKA_path = FUKA_path+'/codes/PythonTools'
pyFUKA_libspath = pyFUKA_path+'/lib/'
sys.path.append(pyFUKA_libspath)
os.environ["HOME_KADATH)"] = FUKA_path
#%env HOME_KADATH=$FUKA_path

# %%
from fuka_plot_tools.setup_utils import get_reader

# %%
interp_var= "gxx"
data_base_path="/home/user/lib/falcon_data"
initial_data_path=pyFUKA_path+'/Example_id/converged_NS_TOTAL_BC.togashi.2.23.-0.4.0.11.dat'
initial_data_path=data_base_path+"/NS_ISO_UNIFORM_ROT.gam2.madm.1.4.omega.0.0.09.dat"
cvar="comega"
res="17"
initial_data_path=data_base_path+f"/NS_ISO_UNIFORM_ROT.gam2.madm.1.4.chi.0.2.0.{res}.dat"
if os.path.isfile(initial_data_path) == True:
    readerISO = get_reader(initial_data_path, ns_iso_uniformrot=True)
else:
    print("{} not found",initial_data_path)

# %%
initial_data_path=data_base_path+"/NS_TOTAL_BC.gam2.madm.1.4.chi.0.0.09.dat"
# initial_data_path=data_base_path+"/NS_TOTAL_BC.gam2.madm.1.4.chi.0.2.0.09.dat"
# print(initial_data_path)
if os.path.isfile(initial_data_path) == True:
    readerNS = get_reader(initial_data_path, ns=True)
else:
    print("{} not found",initial_data_path)

# %%
readerISO.getExporterKeys()

# %%
y1, y2 = [-30, 30]
npts = 512

y_coords = np.linspace(y1, y2, num=npts)

x1, x2 = [-30,30]
x_coords = np.linspace(x1, x2, num=npts)

coords_lst = [[x,y] for y in y_coords for x in x_coords]

# %%
dataISO = readerISO.getFieldValues(cvar, coords_lst, -1)

# %%
data = np.array(dataISO)

# %%
data=data.reshape(len(x_coords),len(x_coords))
# data=np.log(np.abs(data) + 1e-15)

# %%
X, Y = np.meshgrid(x_coords, y_coords)

# %% [markdown]
# 

# %%
norm=None
# norm = colors.Normalize(vmin=-12, vmax = -4, clip=False) 
norm = colors.LogNorm(vmin=1e-15, vmax = 1e-4, clip=False) 
plt.rcParams.update({
'figure.figsize'    : [4, 4],
'image.cmap' : 'inferno'
})
fig,axs = plt.subplots(1,1)
axs=[axs]
cs=axs[0].pcolor(X, Y, np.abs(data)+1e-15, norm=norm)
cbarax,kw = mpl.colorbar.make_axes(axs[0], location='right', orientation='vertical', aspect=30)
cbar = fig.colorbar(
        cs,
        cax=cbarax, ticklocation='right',
        orientation='vertical', extend='both',
#         ticks=np.linspace(norm.vmin, norm.vmax, num=5)
      )


plt.savefig(f"/home/user/lib/fuka_dev/codes/PythonTools/{cvar}-{res}.png", dpi=300, bbox_inches="tight")
plt.show()
quit()
# %%
dataISO = readerISO.getExporterFieldValues(interp_var, coords_lst)
dataISO = np.array(dataISO)
dataISO = dataISO.reshape(npts,npts)

# %% [markdown]
# for v in dataISO:
#     if np.isnan(v.any()):
#         print(v)

# %% [markdown]
# dataNS = readerNS.getExporterFieldValues(interp_var, coords_lst)
# dataNS = np.array(dataNS)
# dataNS = dataNS.reshape(npts,npts)

# %% [markdown]
# X, Y = np.meshgrid(x_coords, y_coords)

# %% [markdown]
# norm = None 
# # norm = colors.LogNorm(vmin=1e-10, vmax = 1e-2, clip=False) 
# norm = colors.Normalize(vmin=1.3, vmax = 2, clip=False) 

# %% [markdown]
# plt.rcParams.update({
# 'figure.figsize'    : [8.0 * 2, 4],
# 'image.cmap' : 'inferno'
# })
# fig,axs = plt.subplots(1,2)
# cs=axs[0].pcolor(X, Y, dataISO, norm=norm)
# axs[1].pcolor(X, Y, dataNS, norm=norm)
# cbarax,kw = mpl.colorbar.make_axes(axs, location='right', orientation='vertical', aspect=30)
# cbar = fig.colorbar(
#         cs,
#         cax=cbarax, ticklocation='right',
#         orientation='vertical', extend='both',
#         ticks=np.linspace(norm.vmin, norm.vmax, num=5)
#       )

# %%
vals = readerISO.getallExporterFieldValues_pointwise([0., 0., 0.])

# %%
def detgij(vals):
    return vals["gxx"] * vals["gyy"] * vals["gzz"] 
+ vals["gxy"] * vals["gyz"] * vals["gxz"] 
+ vals["gxz"] * vals["gxy"] * vals["gyz"] 
- vals["gxz"] * vals["gyy"] * vals["gxz"] 
- vals["gxy"] * vals["gxy"] * vals["gzz"] 
- vals["gxx"] * vals["gyz"] * vals["gyz"]

# %%
detgij(vals)

# %%
for i in np.linspace(0,30, num=256):
    x = 0.
    y = 0.
    z = i
    vals = readerISO.getallExporterFieldValues_pointwise([x,y,z])
    print(detgij(vals), vals["gxx"] , vals["gyy"] , vals["gzz"],
          vals["gxy"] * vals["gxz"] * vals["gyz"])

# %%



