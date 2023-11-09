"""
Designed to process time series data collected on a Zeiss confocal microscope from microfluidic flow experiments.

Data collection: A line scan is selected at one point along the microfluidics channel; the line scan is repeated 
to create a time series experiment. Time series experiments are repeated for a control (in which the fluorescence 
intensity profile should be flat), followed by 2 experiments in which the fluorescence intensity profile is 
expected to be perturbed. The experiments are obtained at different flow rates, to lead to different residence 
times at the same point along the microfluidics channel.

Data storage: Each flow rate has it’s own directory, inside which the 3 time series experiments (control, and two 
others) are saved as CZI files. The control always appears at the top of the direction, alphabetically, by 
appropriate file naming (important for processing).

Using script: Run Cell 0 to set up processing for all flow rates (defining lists: dsets and FR is important for 
comparing across flow rates); only run once at the beginning. Import and check the averaged data across each time 
series in a single directory (flow rate) by running Cell 1, and normalizing with Cell 2. If the data is acceptable, 
use Cell 3 to store the data and flow rate (in dsets and FR, respectively). Repeat Cells 1-3 for each directory. 
Once all of the directories of interest have been imported, checked, normalized and stored (using Cells 1-3), run 
Cell 4 to plot a comparison of the results across directories, which in this case comprises of results at different 
flow rates.
"""



# (0) Set up for processing across flow rates
# run only once at the beginning of processing experiments collected in same sitting
import os
import czifile as cz 
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
import string

#calculate the residence time at given flow rate (in uL/min) and disance from inlets (cm)
def resTime(flowrate,L,w=360,h=50):
    """
    Fuction to calculate residence time of solution in microfluidic channel, 
    given the flow rate the solution is pumped at, and the dimensions
    of the microfluidic channel.

    Parameters
    ----------
    flowrate : float
        flow rate that solution is pumped into device at.
    L : float
        length along channel where images were collected (in cm).
    w : float, optional
        width of channel (in cm). The default is 360.
    h : float, optional
        height of channel (in cm). The default is 50.

    Returns
    -------
    resTime : float
        residence time of the solution in the channel at the point
        where the images are collected.

    """
    Lm = L*0.01 #(cm to m)
    wm = w*1e-6 #(um to m) 
    hm = h*1e-6 #(um to m)
    volume = Lm*wm*hm
    vol_um = volume*1e9 #(m3 to uL)
    FR = flowrate/60 #(uL/min to uL/s)
    resTime = vol_um / FR #in seconds
    return resTime

# pixelsize should be noted from Zeiss software during measurement
pixelsize = input('pixel size (in um)? ')
PS = float(pixelsize)
# pixelsize should be noted from Zeiss software during measurement
L = input('distance from inlets (in cm)? ')
L = float(L)

# defining lists to collect data processed at multiple flow rates (in different folders)
dsets = []
FR = []

#%%
# (1) Inital data import and processing; plot
# opens file explorer to select directory containing data of interest
tk.Tk().withdraw()
path = filedialog.askdirectory()

# process all images in folder, assign them the same name as in the folder
stacknames = ()
for entry in os.scandir(path):
    stack = cz.imread(entry)
    # name each stack with filename (excluding .czi suffix)
    globals()[entry.name[:-4]] = stack 
    stacknames += (entry.name[:-4],) 


# average over each experiment in the time series
# plot the result compared to the raw data
for n,name in enumerate(stacknames):
    allx = np.zeros(len(globals()[name][0,0,0,0,0,:,0]))
    ns = globals()[name].shape[2] # shape[2] is number of cycles in time series
    npix = globals()[name].shape[5] # shape[5] is number of pixels in scan
    for i in range(ns):
        allx = np.add(allx,globals()[name][0,0,i,0,0,:,0])
    globals()[name + '_avg'] = allx/ns # save mean as name_avg


# comparison plot of averaged fluorescence intensity (FI) of every experiment in directory
plt.figure(figsize=(5,3))
for n,name in enumerate(stacknames):
    globals()[name + '_avg'] = pd.DataFrame(globals()[name + '_avg'])
    plt.plot(globals()[name + '_avg'].index*PS,globals()[name + '_avg'],label=name)
plt.xlabel(r'Distance / $\mu$m')
plt.ylabel('FI / a.u.')
plt.xlim(0,PS*512)
plt.legend(frameon=False,prop={'size':9},ncol=2)
#plt.savefig('Example.pdf',bbox_inches='tight')


#%%
# (2) Select region unaffected by microfluidic channel walls, scale and replot
# ACTION REQUIRED: define distance to leave at either side to avoid edge effects
# this may differ for experiments taken on different days/positions
cutL = 30
cutR = 40

control = globals()[stacknames[0] + '_avg'].iloc[cutL:-cutR]
totFI = control.sum()

# normalise FI for every experiment in directory, collect and plot
control = globals()[stacknames[0] + '_avg'].iloc[cutL:-cutR] # define control experiment (flat FI)
totFI = control.sum() # define total FI
ndata = [] # collect normalied dataset for each dataset
plt.figure(figsize=(5,3))
for n,name in enumerate(stacknames):
    data = globals()[name + '_avg'].iloc[cutL:-cutR]
    flattened = data/control # accounts for non-uniformity of laser
    scale = (data/control).sum() / (control/control).sum() # scale so each FI profile has same integration
    ndata_i = flattened/scale
    ndata.append(pd.DataFrame(ndata_i.values, index=(data.index-512/2)*PS)) # adjust index to plot distance from center
    plt.plot((ndata_i.index-512/2)*PS,ndata_i,label=name)
plt.xlabel(r'Distance from center / $\mu$m')
plt.ylabel('NFI / a.u.')
#plt.ylim(0.7,1.2)
plt.legend(frameon=False,prop={'size':9},ncol=2)
#plt.savefig('Example_norm.pdf',bbox_inches='tight')


#%%
# (3) Storing processed data, before repeating processing for other flow rates
# storing result from each flow rate, as calcuated in cell above
dset = input('name of dset? ')
globals()[dset] = ndata 
dsets.append(dset)
# storing flow rate
flrate = input('flow rate? ')
flrate = float(flrate)
FR.append(flrate)

# Need to repeat cells 1-3 for each flow rate before moving to cell 4; do not run cell 0 again

#%%
# (4) Comparing results across flow rates

# PLOT 1
# set appropriate axis titles
title_L = '1 M NaCl'
title_R = r'1 M Na$_2$SO$_4$'
# plot of FI across channel, each flow rate represented by shades of grey
fig, ax = plt.subplots(1,2,figsize=(8,3))
for axs in ax: 
    axs.plot(globals()[dsets[0]][0],'k') # same control on each (all normalised to y=1)
    axs.set_xlabel(r'Distance from center / $\mu$m')
    axs.set_ylabel('NFI / a.u.')
    axs.set_ylim(0.9,1.15)
cmap = plt.get_cmap('Greys')
cols = [cmap(i) for i in np.linspace(0, 1, len(dsets)+3)]
for i,col in zip(range(len(dsets)),cols[3:]):
    ax[0].plot(globals()[dsets[i]][1],color=col)
    ax[1].plot(globals()[dsets[i]][2],color=col)
# ACTION REQUIRED: set appropriate axis titles
ax[0].set_title(title_L)
ax[1].set_title(title_R)
plt.tight_layout()

# PLOT 2
# plot of FI across channel, each flow rate represented by shades of grey (top)
# plot of FI at each residence time, with FI represented on red-to-blue colourmap (bottom)
fig, ax = plt.subplots(2,3,figsize=(8,4.5),gridspec_kw={'width_ratios': [19,19,1]},sharex=('col'),)
ax = ax.flatten()
for axs in [ax[3],ax[4]]: 
    axs.set_ylabel(r'Residence time / s')
    axs.set_xlabel(r'Distance from center / $\mu$m')
    axs.set_xlim(-136,128)
for n,axs in enumerate([ax[0],ax[1]]): 
    axs.text(-0.16,1.05,'('+string.ascii_lowercase[n]+')',transform=axs.transAxes)
    axs.set_ylabel('NFI / a.u.')
    axs.plot(globals()[dsets[0]][0],'k') # can plot same control on each because all normalised to y=1
    axs.set_ylim(0.9,1.14)
    axs.set_xlim(-136,128)
cmap = plt.get_cmap('Greys')
cols = [cmap(i) for i in np.linspace(0, 1, len(dsets)+3)]
cols = cols[3:]
RT = []
for i,col in zip(range(len(dsets)),cols):
    RT.append(resTime(FR[i],L=L))
    ax[0].plot(globals()[dsets[i]][1],color=col)
    ax[1].plot(globals()[dsets[i]][2],color=col)
ax[0].set_title(title_L)
ax[1].set_title(title_R)

# find overall max and min to scale the data appropraitely for colormapping
zmin = []
zmax = []
for i in range(len(dsets)):
    zmin.append(min(min(globals()[dsets[i]][1].iloc[:,0]),min(globals()[dsets[i]][2].iloc[:,0])))
    zmax.append(max(max(globals()[dsets[i]][1].iloc[:,0]),max(globals()[dsets[i]][2].iloc[:,0])))
for axs, n in zip([ax[3],ax[4]],[1,2]):
    for i in range(len(dsets)):
        x = globals()[dsets[i]][n].index
        y = resTime(FR[i],L=L)*np.ones_like(x)
        z = globals()[dsets[i]][n].iloc[:,0]
        # from: https://stackoverflow.com/questions/20165169/change-colour-of-curve-according-to-its-y-value-in-matplotlib
        axs.scatter(x,y, c=cm.bwr((z-min(zmin))/(max(zmax)-min(zmin))), marker=3)

# set up colour bars:
#   for top plots
cmap = mpl.colors.ListedColormap(cols) # reassign colourmap as the section of colours used (excluding the white bits cut off)
cb1 = mpl.colorbar.ColorbarBase(ax[2], cmap=cmap, norm=mpl.colors.BoundaryNorm(boundaries=RT, ncolors=len(cmap.colors)))
cb1.set_label('Residence time / s')
#   for bottom plots
cb2 = mpl.colorbar.ColorbarBase(ax[5], cmap='bwr', norm=mpl.colors.Normalize(vmin=min(zmin),vmax=max(zmax)))
cb2.set_label(r'NFI / a.u.')

plt.tight_layout()
#plt.savefig('500mM_salt_cbar.pdf',bbox_inches='tight')

