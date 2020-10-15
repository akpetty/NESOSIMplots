""" plot_nesosim_oib_comp_multiyears.py
	
	Script to compare/plot the already binned OIB data to the coincident 
	NESOSIM v1 data for a given year then all years of that range
	
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: Binned OIB snow depths, NESOSIM data
	Output: Scatter plots and map sanity checks

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		10/10/2020: Version 1

"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import utils as ut
import os
import pyproj
import cartopy.crs as ccrs

from config import forcing_save_path
from config import model_save_path
from config import figure_path
nesosim_og_path = '/Users/aapetty/Data/NESOSIMv1/'

def mask_data_2d_1d(nesosim_data, oib_data):
	# mask data using common filters

	maskDay=np.zeros((nesosim_data.shape[0], nesosim_data.shape[1]))
	maskDay[np.where(np.isnan(nesosim_data))]=1
	maskDay[np.where(np.isnan(oib_data))]=1

	maskDay[np.where(oib_data<=4)]=1
	maskDay[np.where(nesosim_data<=4)]=1

	maskDay[np.where(oib_data>80)]=1
	maskDay[np.where(nesosim_data>80)]=1

	# Apply mask
	nesosim_data= nesosim_data[maskDay<1]
	oib_data= oib_data[maskDay<1]

	return nesosim_data, oib_data


dx_oib=100000
snowType='MEDIAN'
extraStr='v11'
outStr='rc_mask'

mask_region=True
mask_coast=True

dxStr_oib=str(int(dx_oib/1000))+'km'
xptsG, yptsG, latG, lonG, proj = ut.create_grid(dxRes=dx_oib)

month1=7 # 8=September
day1=14
month2=3 # 4=May
day2=29

years=np.arange(2011, 2015+1)

nesosim_1d_all=[]
oib_1d_all=[]

for x in range(len(years)):
	oib_year=years[x]
	oib_path= forcing_save_path+dxStr_oib+'/OIB/'+str(oib_year)+'/'
	files = glob(oib_path+'*'+dxStr_oib+snowType+extraStr)

	dates_oib_year=[file.split('/')[-1][0:8] for file in files]
	print(dates_oib_year)

	_, _, _, dateOut=ut.getDays(oib_year-1, month1, day1, oib_year, month2, day2)


	OIBdaysAfterSepT= ut.getOIBbudgetDays(dates_oib_year, oib_year-1, month1, day1)

	nesosim_1d=[]
	oib_1d=[]
	nesosim_grids=[]
	oib_grids=[]

	for day_oib in range(np.size(dates_oib_year)):
		
		snowDepthMv1, latv1, lonv1 = ut.get_nesosim_day(nesosim_og_path, oib_year, OIBdaysAfterSepT[day_oib], converttocm=1)
		xptsv1, yptsv1 = proj(lonv1, latv1)

		# Linearly interpoalte the model to the 100 km nesosim grid
		nesosim_day_grid = griddata((xptsv1.flatten(), yptsv1.flatten()), snowDepthMv1.flatten(), (xptsG, yptsG), method='linear')

		oib_day_grid=np.load(files[day_oib], allow_pickle=True)
		oib_day_grid=oib_day_grid*100.

		if mask_region:
			oib_day_grid= ma.masked_where(region_maskG>8.5, oib_day_grid)
			nesosim_day_grid= ma.masked_where(region_maskG>8.5, nesosim_day_grid)

		if mask_coast:
			oib_day_grid= ma.masked_where(coast_maskG<50, oib_day_grid)
			nesosim_day_grid= ma.masked_where(coast_maskG<50, nesosim_day_grid)

		nesosim_day_1d, oib_day_1d = mask_data_2d_1d(nesosim_day_grid, oib_day_grid)
		
		nesosim_grids.append(nesosim_day_grid)
		oib_grids.append(oib_day_grid)

		nesosim_1d.extend(nesosim_day_1d)
		oib_1d.extend(oib_day_1d)

	nesosim_grid=np.nanmean(nesosim_grids, axis=0)
	oib_grid=np.nanmean(oib_grids, axis=0)

	nesosim_1d_all.extend(nesosim_1d)
	oib_1d_all.extend(oib_1d)


	plot_map=True
	if (plot_map):
		fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3.7))

		ax=axs.flatten()[0]
		plt.sca(ax)
		plt.imshow(nesosim_grid, origin='lower', vmin=0, vmax=60)

		ax=axs.flatten()[1]
		plt.sca(ax)
		plt.imshow(oib_grid, origin='lower', vmin=0, vmax=60)

		plt.savefig(figure_path+'/OIB/comps/'+str(oib_year)+'OIBcorrelations_nesosimv1_'+dxStr_oib+snowType+outStr+'_map.png', dpi=300)


	fig = plt.figure(figsize=(4, 3))
	ax=plt.gca()
	plt.sca(ax)

	im1 = plt.scatter(nesosim_1d,oib_1d, color='0.4',s=6, marker='x', alpha=0.7)

	plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

	trend, sig, r_a, intercept = ut.correlateVars(nesosim_1d,oib_1d)
	r_str = '%.2f' % r_a
	rmse=np.sqrt(np.mean((np.array(nesosim_1d)-np.array(oib_1d))**2))
	rmsStr='%.0f' % rmse

	ax.annotate(str(oib_year), xy=(0.02, 1.01), 
			xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	ax.annotate('r: '+r_str+'  RMSE: '+rmsStr+' cm', xy=(0.02, 0.92), 
			xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	plt.xlim(0, 80)
	plt.ylim(0, 80)

	plt.ylabel('OIB snow depth '+snowType+' (cm)')
	plt.xlabel('NESOSIMv1 snow depth (cm)') 

	plt.tight_layout()
	plt.savefig(figure_path+'/OIB/comps/'+str(oib_year)+'OIBcorrelations_nesosimv1_'+dxStr_oib+snowType+outStr+'.pdf', dpi=300)
	#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
	plt.close(fig)


fig = plt.figure(figsize=(4, 3))
ax=plt.gca()
plt.sca(ax)

im1 = plt.scatter(nesosim_1d_all,oib_1d_all, color='0.4',s=6, marker='x', alpha=0.7)

plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

trend, sig, r_a, intercept = ut.correlateVars(nesosim_1d_all,oib_1d_all)
r_str = '%.2f' % r_a
rmse=np.sqrt(np.mean((np.array(nesosim_1d_all)-np.array(oib_1d_all))**2))
rmsStr='%.0f' % rmse

ax.annotate(str(years[0])+'-'+str(years[-1]), xy=(0.02, 1.01), 
		xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

ax.annotate('r: '+r_str+'  RMSE: '+rmsStr+' cm', xy=(0.02, 0.92), 
		xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

plt.xlim(0, 80)
plt.ylim(0, 80)

plt.ylabel('OIB snow depth '+snowType+' (cm)')
plt.xlabel('NESOSIMv1 snow depth (cm)') 

plt.tight_layout()
plt.savefig(figure_path+'/OIB/comps/'+str(years[0])+'-'+str(years[-1])+'OIBcorrelations_nesosimv1_'+dxStr_oib+snowType+outStr+'.pdf', dpi=300)
#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
plt.close(fig)



