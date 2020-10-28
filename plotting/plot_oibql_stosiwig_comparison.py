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

def mask_data_2d_1d(data1, data2):
	""" mask data using common filters"""


	maskDay=np.zeros((data1.shape[0], data1.shape[1]))
	maskDay[np.where(~np.isfinite(data1))]=1
	maskDay[np.where(~np.isfinite(data2))]=1

	# Mask less than 4 cm
	maskDay[np.where(data2<=3)]=1
	maskDay[np.where(data1<=3)]=1

	# Mask greater than 80 cm
	maskDay[np.where(data2>80)]=1
	maskDay[np.where(data1>80)]=1

	# Apply mask
	data1= data1[maskDay<1]
	data2= data2[maskDay<1]

	return data1, data2

dx_oib=100000
extraStr='v11'
outStr='rc_mask'
snowType1=''
snowType2='MEDIAN'
mask_region=False
mask_coast=False

dxStr_oib=str(int(dx_oib/1000))+'km'
xptsG, yptsG, latG, lonG, proj = ut.create_grid(dxRes=dx_oib)

region_mask, xptsI, yptsI = ut.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

xptsc, yptsc, coast_mask = ut.get_coastal_mask(anc_data_path, proj)
coast_maskG = griddata((xptsc.flatten(), yptsc.flatten()), coast_mask.flatten(), (xptsG, yptsG), method='nearest')

month1=8 # 8=September
day1=0
month2=3 # 4=May
day2=29

years=np.arange(2010, 2015+1)

oib_1d_all1=[]
oib_1d_all2=[]

for x in range(len(years)):
	oib_year=years[x]
	oib_path= forcing_save_path+dxStr_oib+'/OIB/'+str(oib_year)+'/'
	files1 = glob(oib_path+'*'+dxStr_oib+snowType1+extraStr)
	files2 = glob(oib_path+'*'+dxStr_oib+snowType2+extraStr)

	dates_oib_year1=[file.split('/')[-1][0:8] for file in files1]
	dates_oib_year2=[file.split('/')[-1][0:8] for file in files2]
	print(dates_oib_year1)
	print(dates_oib_year2)

	# just get dates where they match
	dates_oib_year=[date for date in dates_oib_year1 if date in dates_oib_year2]
	print(dates_oib_year)

	_, _, _, dateOut=ut.getDays(oib_year-1, month1, day1, oib_year, month2, day2)

	OIBdaysAfterSepT= ut.getOIBbudgetDays(dates_oib_year, oib_year-1, month1, day1)

	oib_1d1=[]
	oib_1d2=[]
	oib_grids1=[]
	oib_grids2=[]

	for date in dates_oib_year:
		
		file1 = next((s for s in files1 if date in s), None)
		print(file1)
		file2 = next((s for s in files2 if date in s), None)
		print(file2)

		oib_day_grid1=np.load(file1, allow_pickle=True)
		oib_day_grid1= ma.masked_where(np.isnan(oib_day_grid1), oib_day_grid1)

		oib_day_grid2=np.load(file2, allow_pickle=True)
		oib_day_grid2= ma.masked_where(np.isnan(oib_day_grid2), oib_day_grid2)

		if mask_region:
			oib_day_grid1= ma.masked_where(region_maskG>8.5, oib_day_grid1)
			oib_day_grid2= ma.masked_where(region_maskG>8.5, oib_day_grid2)


		if mask_coast:
			oib_day_grid1 = ma.masked_where(coast_maskG<50, oib_day_grid1)
			oib_day_grid2 = ma.masked_where(coast_maskG<50, oib_day_grid2)

		oib_day_grid1=oib_day_grid1*100.
		oib_day_grid2=oib_day_grid2*100.

		oib_day_1d1, oib_day_1d2 = mask_data_2d_1d(oib_day_grid1, oib_day_grid2)
		print(np.size(oib_day_1d1))

		oib_grids1.append(oib_day_grid1)
		oib_grids2.append(oib_day_grid2)

		oib_1d1.extend(oib_day_1d1)
		oib_1d2.extend(oib_day_1d2)

		plot_map=True
		if (plot_map):
			fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 3.5))

			ax=axs.flatten()[0]
			plt.sca(ax)
			plt.imshow(oib_day_grid1, origin='lower', vmin=0, vmax=60)

			ax=axs.flatten()[1]
			plt.sca(ax)
			plt.imshow(oib_day_grid2, origin='lower', vmin=0, vmax=60)

			ax=axs.flatten()[2]
			plt.sca(ax)
			plt.imshow(oib_day_grid2-oib_day_grid1, origin='lower', vmin=-30, vmax=30)


			plt.savefig(figure_path+'/OIB/ql_comps/'+str(oib_year)+'OIBcorrelations'+dxStr_oib+extraStr+snowType1+snowType2+date+'_map.png', dpi=300)



	oib_grid1=ma.mean(oib_grids1, axis=0)
	oib_grid2=ma.mean(oib_grids2, axis=0)

	oib_1d_all1.extend(oib_1d1)
	oib_1d_all2.extend(oib_1d2)

	plot_map=True
	if (plot_map):
		fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 3.5))

		ax=axs.flatten()[0]
		plt.sca(ax)
		plt.imshow(oib_grid1, origin='lower', vmin=0, vmax=60)

		ax=axs.flatten()[1]
		plt.sca(ax)
		plt.imshow(oib_grid2, origin='lower', vmin=0, vmax=60)

		ax=axs.flatten()[2]
		plt.sca(ax)
		plt.imshow(oib_grid2-oib_grid1, origin='lower', vmin=-30, vmax=30)


		plt.savefig(figure_path+'/OIB/ql_comps/'+str(oib_year)+'OIBcorrelations'+dxStr_oib+extraStr+snowType1+snowType2+'_map.png', dpi=300)


	fig = plt.figure(figsize=(4, 3))
	ax=plt.gca()
	plt.sca(ax)

	im1 = plt.scatter(oib_1d1,oib_1d2, color='0.4',s=6, marker='x', alpha=0.7)

	plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

	trend, sig, r_a, intercept = ut.correlateVars(oib_1d1,oib_1d2)
	#r_str = '%.2f' % r_a
	rmse=np.sqrt(np.mean((np.array(oib_1d1)-np.array(oib_1d2))**2))
	rmsStr='%.0f' % rmse

	rStr = '%.2f' % r_a

	merr=np.mean(np.array(oib_1d1)-np.array(oib_1d2))
	merrStr='%.1f' % merr	

	std=np.std(np.array(oib_1d1)-merr-np.array(oib_1d2))
	stdStr='%.1f' % std

	num_grids=np.size(oib_1d1)

	ax.annotate(str(oib_year), xy=(0.02, 1.01), 
			xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	#ax.annotate('r: '+r_str+'  RMSE: '+rmsStr+' cm', xy=(0.02, 0.92), 
	#		xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	ax.annotate('r: '+rStr+'\nRMSE: '+rmsStr+' cm'+'\nMean bias: '+merrStr+' cm'+'\nSD: '+stdStr+' cm'+'\nN: '+str(num_grids), xy=(0.02, 0.98), 
			xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')

	plt.xlim(0, 80)
	plt.ylim(0, 80)

	plt.ylabel('OIB ('+snowType2+') snow depth (cm)')
	plt.xlabel('OIB QL snow depth (cm)') 

	plt.tight_layout()
	plt.savefig(figure_path+'/OIB/ql_comps/'+str(oib_year)+'OIBcorrelations'+dxStr_oib+extraStr+snowType1+snowType2+'.pdf', dpi=300)
	#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
	plt.close(fig)


fig = plt.figure(figsize=(4, 3))
ax=plt.gca()
plt.sca(ax)

im1 = plt.scatter(oib_1d_all1,oib_1d_all2, color='0.4',s=6, marker='x', alpha=0.7)

plt.plot(np.arange(0, 100, 0.1), np.arange(0, 100, 0.1), 'k', ls='--')

trend, sig, r_a, intercept = ut.correlateVars(oib_1d_all1,oib_1d_all2)

rStr = '%.2f' % r_a

rmse=np.sqrt(np.mean((np.array(oib_1d_all1)-np.array(oib_1d_all2))**2))
rmsStr='%.0f' % rmse

merr=np.mean(np.array(oib_1d_all1)-np.array(oib_1d_all2))
merrStr='%.1f' % merr	

std=np.std(np.array(oib_1d_all1)-merr-np.array(oib_1d_all2))
stdStr='%.1f' % std

num_grids=np.size(oib_1d_all1)

ax.annotate(str(years[0])+'-'+str(years[-1]), xy=(0.02, 1.01), 
		xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

ax.annotate('r: '+rStr+'\nRMSE: '+rmsStr+' cm'+'\nMean bias: '+merrStr+' cm'+'\nSD: '+stdStr+' cm'+'\nN: '+str(num_grids), xy=(0.02, 0.98), 
		xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')


plt.xlim(0, 80)
plt.ylim(0, 80)

plt.ylabel('OIB ('+snowType2+') snow depth (cm)')
plt.xlabel('OIB QL snow depth (cm)') 

plt.tight_layout()
plt.savefig(figure_path+'/OIB/ql_comps/'+str(years[0])+'-'+str(years[-1])+'OIBcorrelations'+dxStr_oib+extraStr+snowType1+snowType2+'.pdf', dpi=300)
#savefig(figpath+'/seasonalSnowDensityComp4'+folderStr+'.png', dpi=300)
plt.close(fig)

