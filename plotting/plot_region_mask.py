
""" plot_region_mask.py
	
	Python dependencies:
		See below for the relevant module imports
		Also reads in some functions in utils.py

"""


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import utils as cF
import os
import pyproj
import cartopy.crs as ccrs

from config import reanalysis_raw_path
from config import forcing_save_path
from config import figure_path

anc_data_path='../../anc_data/'

dx=50000
xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
print(xptsG)
print(yptsG)

dxStr=str(int(dx/1000))+'km'
print(dxStr)


region_mask, xptsI, yptsI = cF.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

cF.plot_gridded_cartopy(lonG, latG, region_maskG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/region_mask', varStr='Region ', units_lab=r'#', minval=0, maxval=12, cmap_1=plt.cm.viridis)
		