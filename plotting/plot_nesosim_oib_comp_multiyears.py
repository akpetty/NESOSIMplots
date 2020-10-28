""" plot_nesosim_oib_comp_multiyears.py
	
	Script to compare/plot the already binned OIB data to the coincident 
	NESOSIM v11 (new version) data for a given year then all years of that range
	
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

from config import model_save_path, figure_path, forcing_save_path, anc_data_path

def mask_data_2d_1d(nesosim_data, oib_data):
	""" mask data using common filters"""

	maskDay=np.zeros((nesosim_data.shape[0], nesosim_data.shape[1]))
	maskDay[np.where(np.isnan(nesosim_data))]=1
	maskDay[np.where(np.isnan(oib_data))]=1

	maskDay[np.where(oib_data<=3)]=1
	maskDay[np.where(nesosim_data<=3)]=1

	maskDay[np.where(oib_data>80)]=1
	maskDay[np.where(nesosim_data>80)]=1

	# Apply mask
	nesosim_data= nesosim_data[maskDay<1]
	oib_data= oib_data[maskDay<1]

	return nesosim_data, oib_data

dx_oib=100000
extraStr='v11'
outStr='rc_mask'
snowType=''
mask_region=True
mask_coast=True

dxStr_oib=str(int(dx_oib/1000))+'km'
xptsG, yptsG, latG, lonG, proj = ut.create_grid(dxRes=dx_oib)

dx_m=100000
dxStr_m=str(int(dx_m/1000))+'km'
if (dx_m<dx_oib):
	xptsGm, yptsGm, latGm, lonGm, projm = ut.create_grid(dxRes=dx_m)

region_mask, xptsI, yptsI = ut.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

xptsc, yptsc, coast_mask = ut.get_coastal_mask(anc_data_path, proj)
coast_maskG = griddata((xptsc.flatten(), yptsc.flatten()), coast_mask.flatten(), (xptsG, yptsG), method='nearest')


month1=8 # 8=September
day1=0
month2=3 # 4=May
day2=29

years=np.arange(2010, 2015+1)

nesosim_1d_all=[]
oib_1d_all=[]

for x in range(len(years)):
	oib_year=years[x]
	oib_path= forcing_save_path+dxStr_oib+'/OIB/'+str(oib_year)+'/'
	files = glob(oib_path+'*'+dxStr_oib+snowType+extraStr)

	dates_oib_year=[file.split('/')[-1][0:8] for file in files]
	print(dates_oib_year)

	# get rid of dates beyong Apr 30

	dates_oib_year=[date for date in dates_oib_year if date[4:6] not in '05']
	print(dates_oib_year)

	_, _, _, dateOut=ut.getDays(oib_year-1, month1, day1, oib_year, month2, day2)

	OIBdaysAfterSepT= ut.getOIBbudgetDays(dates_oib_year, oib_year-1, month1, day1)


	reanalysis='ERA5'
	folderStr=reanalysis+'sfERA5windsOSISAFdriftsCDRsicrhovariable_IC2_DYN1_WP1_LL1_AL1_WPF5.8e-07_WPT5_LLF1.45e-07-'+dxStr_m+'v11oct28'
	totalOutStr=''+folderStr+'-'+dateOut


	nesosim_1d=[]
	oib_1d=[]
	nesosim_grids=[]
	oib_grids=[]

	for day_oib in range(np.size(dates_oib_year)):
		nesosim_day_grid=ut.get_budgets2layers_day(['snowDepthTotalConc'], model_save_path+dxStr_m+'/', folderStr, OIBdaysAfterSepT[day_oib], totalOutStr)
		oib_day_grid=np.load(files[day_oib], allow_pickle=True)
		oib_day_grid= ma.masked_where(np.isnan(oib_day_grid), oib_day_grid)

		if mask_region:
			oib_day_grid= ma.masked_where(region_maskG>8.5, oib_day_grid)
			nesosim_day_grid= ma.masked_where(region_maskG>8.5, nesosim_day_grid)

		if mask_coast:
			oib_day_grid= ma.masked_where(coast_maskG<50, oib_day_grid)
			nesosim_day_grid= ma.masked_where(coast_maskG<50, nesosim_day_grid)

		nesosim_day_grid=nesosim_day_grid*100.
		oib_day_grid=oib_day_grid*100.

		if (dx_m<dx_oib):
			# Linearly interpoalte the model to the 100 km nesosim grid
			nesosim_day_grid = griddata((xptsGm.flatten(), yptsGm.flatten()), nesosim_day_grid.flatten(), (xptsG, yptsG), method='linear')

		nesosim_day_1d, oib_day_1d = mask_data_2d_1d(nesosim_day_grid, oib_day_grid)
		
		nesosim_grids.append(nesosim_day_grid)
		oib_grids.append(oib_day_grid)

		nesosim_1d.extend(nesosim_day_1d)
		oib_1d.extend(oib_day_1d)

	nesosim_grid=ma.mean(nesosim_grids, axis=0)
	oib_grid=ma.mean(oib_grids, axis=0)

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

		plt.savefig(figure_path+'/OIB/model_comps/'+str(oib_year)+'OIBcorrelations'+folderStr+dxStr_m+dxStr_oib+extraStr+snowType+outStr+'_map.png', dpi=300)


	fig = plt.figure(figsize=(4, 3))
	ax=plt.gca()
	plt.sca(ax)

	im1 = plt.scatter(nesosim_1d,oib_1d, color='0.4',s=6, marker='x', alpha=0.7)

	plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

	trend, sig, r_a, intercept = ut.correlateVars(nesosim_1d,oib_1d)
	rStr = '%.2f' % r_a

	rmse=np.sqrt(np.mean((np.array(nesosim_1d)-np.array(oib_1d))**2))
	rmsStr='%.0f' % rmse

	merr=np.mean(np.array(nesosim_1d)-np.array(oib_1d))
	merrStr='%.1f' % merr	

	std=np.std(np.array(nesosim_1d)-merr-np.array(oib_1d))
	stdStr='%.1f' % std

	num_grids=np.size(nesosim_1d)

	ax.annotate(str(oib_year), xy=(0.02, 1.01), 
			xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	ax.annotate('r: '+rStr+'\nRMSE: '+rmsStr+' cm'+'\nMean bias: '+merrStr+' cm'+'\nSD: '+stdStr+' cm'+'\nN: '+str(num_grids), xy=(0.02, 0.98), 
			xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')

	plt.xlim(0, 80)
	plt.ylim(0, 80)

	plt.ylabel('OIB snow depth '+snowType+' (cm)')
	plt.xlabel('NESOSIM v11 ('+reanalysis+'-SF) snow depth (cm)') 

	plt.tight_layout()
	plt.savefig(figure_path+'/OIB/model_comps/'+str(oib_year)+'OIBcorrelations'+folderStr+dxStr_m+dxStr_oib+extraStr+snowType+outStr+'.pdf', dpi=300)
	#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
	plt.close(fig)


fig = plt.figure(figsize=(4, 3))
ax=plt.gca()
plt.sca(ax)

im1 = plt.scatter(nesosim_1d_all,oib_1d_all, color='0.4',s=6, marker='x', alpha=0.7)

plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

trend, sig, r_a, intercept = ut.correlateVars(nesosim_1d_all,oib_1d_all)
rStr = '%.2f' % r_a

rmse=np.sqrt(np.mean((np.array(nesosim_1d_all)-np.array(oib_1d_all))**2))
rmsStr='%.0f' % rmse

merr=np.mean(np.array(nesosim_1d_all)-np.array(oib_1d_all))
merrStr='%.1f' % merr	

std=np.std(np.array(nesosim_1d_all)-merr-np.array(oib_1d_all))
stdStr='%.1f' % std

num_grids=np.size(nesosim_1d_all)

ax.annotate(str(years[0])+'-'+str(years[-1]), xy=(0.02, 1.01), 
		xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

ax.annotate('r: '+rStr+'\nRMSE: '+rmsStr+' cm'+'\nMean bias: '+merrStr+' cm'+'\nSD: '+stdStr+' cm'+'\nN: '+str(num_grids), xy=(0.02, 0.98), 
			xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')

plt.xlim(0, 80)
plt.ylim(0, 80)

plt.ylabel('OIB snow depth '+snowType+' (cm)')
plt.xlabel('NESOSIM v1.1 snow depth (cm)') 

plt.tight_layout()
plt.savefig(figure_path+'/OIB/model_comps/'+str(years[0])+'-'+str(years[-1])+'OIBcorrelations'+folderStr+dxStr_m+dxStr_oib+extraStr+snowType+outStr+'.pdf', dpi=300)
#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
plt.close(fig)

