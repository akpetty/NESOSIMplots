""" compOIB.py
	
	Script to compare/plot the already binned OIB data to the coincident NESOSIM data
	Model written by Alek Petty (05/22/2019)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: Binned OIB snow depths, NESOSIM data
	Output: Scatter plots

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		05/22/2019: Version 1

"""

import matplotlib
matplotlib.use("AGG")

from glob import glob
import matplotlib.pyplot as plt
import matplotlib.colorbar as mcbar
from scipy.interpolate import griddata
import xarray as xr
from scipy import stats
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import pyproj
import numpy.ma as ma
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import utils as ut


from config import forcing_save_path
from config import model_save_path
from config import figure_path

nesosim_og_path = '/Users/aapetty/Data/NESOSIMv1/'

dx_m=100000
xptsG, yptsG, latG, lonG, projm = ut.create_grid(dxRes=dx_m)
dxStr=str(int(dx_m/1000))+'km'

# STart years (not end OIB years like the other scripts)
years=np.arange(2010, 2014+1)
snowDepthMv1_100kms=[]
snowDepthMv11s=[]

for x in range(len(years)):

	# NESOSIM dates
	month1=8 # 8=September
	day1=0
	month2=3 # 4=May
	day2=29
	yearS=years[x]

	_, _, _, dateOut=ut.getDays(yearS, month1, day1, yearS+1, month2, day2)

	reanalysis='ERA5'
	folderStr=reanalysis+'sfERAIwindsNSIDCdriftsCDRsicrhovariable_IC2_DYN1_WP1_LL0_AL1_WPF5.8e-07_WPT5_LLF2.9e-07-100kmv11oct15'
	totalOutStr=''+folderStr+'-'+dateOut

	# take the last day of both!
	snowDepthMv11=ut.get_budgets2layers_day(['snowDepthTotalConc'], model_save_path+dxStr+'/', folderStr, -1, totalOutStr, converttocm=1)
	snowDepthMv11[~np.isfinite(snowDepthMv11)]=np.nan

	snowDepthMv1, latv1, lonv1 = ut.get_nesosim_day(nesosim_og_path, yearS+1, -1, converttocm=1)
	snowDepthMv1[~np.isfinite(snowDepthMv1)]=np.nan

	xptsv1, yptsv1 = projm(lonv1, latv1)

	# Linearly interpoalte the model to the 100 km nesosim grid
	snowDepthMv1_100km = griddata((xptsv1.flatten(), yptsv1.flatten()), snowDepthMv1.flatten(), (xptsG, yptsG), method='linear')


	snowDepthMv1_100kms.append(snowDepthMv1_100km)
	snowDepthMv11s.append(snowDepthMv11)

	minval=0
	maxval=60

	fig=plt.figure(figsize=(7, 3))
	ax1 = plt.subplot(1, 3, 1, projection=ccrs.NorthPolarStereo(central_longitude=-45))
	plt.subplots_adjust(bottom=0.13, left=0.02, right=0.98, wspace=0.05)

	cs=ax1.pcolormesh(lonG, latG, snowDepthMv1_100km, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

	ax1.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
	ax1.gridlines()
	ax1.coastlines(linewidth=0.25) 
	ax1.set_title('v1 ('+dateOut[-8:]+')')
	cax,kw = mcbar.make_axes(ax1,location='bottom',pad=0.05,shrink=0.7)
	cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
	cb.set_label('snow depth (cm)',size=8)

	#p1.cbar.set_label('CTL run monthly mean skin temperature (C)')

	ax2 = plt.subplot(1, 3, 2, projection=ccrs.NorthPolarStereo(central_longitude=-45))
	cs2=ax2.pcolormesh(lonG, latG, snowDepthMv11, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

	ax2.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
	ax2.gridlines()
	ax2.coastlines(linewidth=0.25)
	ax2.set_title('v2 ('+dateOut[-8:]+')')
	cax2,kw = mcbar.make_axes(ax2,location='bottom',pad=0.05,shrink=0.7)
	cb2=fig.colorbar(cs2,cax=cax2,extend='both',**kw)
	cb2.set_label('snow depth (cm)',size=8)
	#p2.cbar.set_label('2 day relaxation monthly mean skin temperature (C)' )

	ax3 = plt.subplot(1, 3, 3, projection=ccrs.NorthPolarStereo(central_longitude=-45))

	cs3=ax3.pcolormesh(lonG, latG, snowDepthMv11-snowDepthMv1_100km, vmin=-20, vmax=20, cmap=plt.cm.RdBu_r, transform=ccrs.PlateCarree(), zorder=2)

	ax3.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
	ax3.gridlines()
	ax3.coastlines(linewidth=0.25)   
	ax3.set_title('difference')
	cax3,kw = mcbar.make_axes(ax3,location='bottom',pad=0.05,shrink=0.7)
	cb3=fig.colorbar(cs3,cax=cax3,extend='both',**kw)
	cb3.set_label('difference(cm)',size=8)

	plt.savefig(figure_path+'/model_comps/v1_v11'+dateOut[-8:]+folderStr+'.png', dpi=300)

snowDepthMv1_100km_mean=np.nanmean(snowDepthMv1_100kms, axis=0)
snowDepthMv11_mean=np.nanmean(snowDepthMv11s, axis=0)

fig=plt.figure(figsize=(7, 3))
ax1 = plt.subplot(1, 3, 1, projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.subplots_adjust(bottom=0.13, left=0.02, right=0.98, wspace=0.05)

cs=ax1.pcolormesh(lonG, latG, snowDepthMv1_100km_mean, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

ax1.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
ax1.gridlines()
ax1.coastlines(linewidth=0.25) 
ax1.set_title('v1 ('+dateOut[-8:-4]+str(years[0])+'-'+str(years[-1])+')')
cax,kw = mcbar.make_axes(ax1,location='bottom',pad=0.05,shrink=0.7)
cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
cb.set_label('snow depth (cm)',size=8)

#p1.cbar.set_label('CTL run monthly mean skin temperature (C)')

ax2 = plt.subplot(1, 3, 2, projection=ccrs.NorthPolarStereo(central_longitude=-45))
cs2=ax2.pcolormesh(lonG, latG, snowDepthMv11_mean, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

ax2.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
ax2.gridlines()
ax2.coastlines(linewidth=0.25)
ax2.set_title('v2 ('+dateOut[-8:-4]+')')
cax2,kw = mcbar.make_axes(ax2,location='bottom',pad=0.05,shrink=0.7)
cb2=fig.colorbar(cs2,cax=cax2,extend='both',**kw)
cb2.set_label('snow depth (cm)',size=8)
#p2.cbar.set_label('2 day relaxation monthly mean skin temperature (C)' )

ax3 = plt.subplot(1, 3, 3, projection=ccrs.NorthPolarStereo(central_longitude=-45))

cs3=ax3.pcolormesh(lonG, latG, snowDepthMv11_mean-snowDepthMv1_100km_mean, vmin=-20, vmax=20, cmap=plt.cm.RdBu_r, transform=ccrs.PlateCarree(), zorder=2)

ax3.set_extent([-179,179, 45, 90], ccrs.PlateCarree())
ax3.gridlines()
ax3.coastlines(linewidth=0.25)   
ax3.set_title('difference')
cax3,kw = mcbar.make_axes(ax3,location='bottom',pad=0.05,shrink=0.7)
cb3=fig.colorbar(cs3,cax=cax3,extend='both',**kw)
cb3.set_label('difference(cm)',size=8)

plt.savefig(figure_path+'/model_comps/v1_v11'+dateOut[-8:-4]+str(years[0])+'-'+str(years[-1])+folderStr+'.png', dpi=300)



