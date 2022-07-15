import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

dataFile = 'data/fMRI/group/featquery.csv'
dataAll = pd.read_csv(dataFile)

selectivities = {'body': ['eba', 'fba'],
				 'face': ['ofa', 'ffa', 'psts'],
				 'house': ['opa', 'ppa', 'rsc'],
				 'spiky': ['nml_lat_1', 'nml_lat_2', 'nml_med_1', 'nml_med_2']}

shapes = ['o','^','s']

# add columns for region, hemisphere and voxel counts
regions, hemispheres, nVoxels, selectivity = [[], [], [], []]
for row in dataAll.index:
	if dataAll['roi'][row].split(sep='_')[0] != 'nml':
		regions.append(dataAll['roi'][row].split(sep='_')[0])
		hemispheres.append(dataAll['roi'][row].split(sep='_')[1])
		nVoxels.append(int(dataAll['roi'][row].split(sep='_')[2]))
	else:
		if dataAll['roi'][row].endswith('noCatVoxels'):
			regions.append(dataAll['roi'][row][:9]+ '_noCatVoxels')
		else:
			regions.append(dataAll['roi'][row][:9])
		hemispheres.append(dataAll['roi'][row].split(sep='_')[3])
		nVoxels.append(int(dataAll['roi'][row].split(sep='_')[4]))
	if regions[-1] in ['eba', 'fba']:
		selectivity.append('body')
	elif regions[-1] in ['ofa', 'ffa', 'psts']:
		selectivity.append('face')
	elif regions[-1] in ['opa', 'ppa', 'rsc']:
		selectivity.append('house')
	elif regions[-1].startswith('nml'):
		selectivity.append('spiky')

dataAll['region'] = regions
dataAll['hemi'] = hemispheres
dataAll['nVoxels'] = nVoxels
dataAll['selectivity'] = selectivity

for topup in ['topUp']:#, 'noTopUp']:
	for b0 in ['noB0']: #, 'b0']:
		for HRFmodel in ['doubleGamma']: #, 'singleGamma']:
			for voxelSize in [64]: #[16,64,256]:


				# 2D scatter plot on animacy aspect ratio space
				## categoryselective regions
				preferences = {}

				resDir = os.path.join('analysis/results/fMRI/ROIsOnPCspace', topup, b0, HRFmodel, 'allSubjects')
				os.makedirs(resDir, exist_ok=True)
				fig, ax = plt.subplots(figsize=(4,4))

				for selectivity, colour in zip(['body','face','house'],['chartreuse','deepskyblue','red']):
					preferences[selectivity] = {}
					for r, region, in enumerate(selectivities[selectivity]):
						shape = shapes[r]
						preferences[selectivity][region] = {}
						for component, contrast in zip(['animacy', 'aspectRatio'], ['animate>inanimate', 'spiky>stubby']):
							preferences[selectivity][region][component] = []
							for subject in experiment['scanInfo'].keys():
								theseData = dataAll[(dataAll['topup'] == topup) &
													(dataAll['b0'] == b0) &
													(dataAll['HRFmodel'] == HRFmodel) &
													(dataAll['nVoxels'] == voxelSize) &
													(dataAll['region'] == region) &
													(dataAll['condition'] == contrast) &
													(dataAll['subject'] == subject)]
								if len(theseData) > 0:
									if component == 'animacy':
										preferences[selectivity][region][component].append(-np.mean(theseData['PSC'])) #  axis inverted
									elif component == 'aspectRatio':
										preferences[selectivity][region][component].append(np.mean(theseData['PSC']))
								else:
									preferences[selectivity][region][component].append(np.nan)

						ax.scatter(preferences[selectivity][region]['aspectRatio'], preferences[selectivity][region]['animacy'], s=30, color=colour, marker=shape)
				for selectivity, colour in zip(['body','face','house'],['chartreuse','deepskyblue','red']):
					for r, region, in enumerate(selectivities[selectivity]):
						shape = shapes[r]
						ax.scatter(np.nanmean(preferences[selectivity][region]['aspectRatio']), np.nanmean(preferences[selectivity][region]['animacy']), s=100,
								   edgecolors='black', marker=shape, color=colour, label=region)
				ax.set_title(f'category selective ROIs in object space')
				plt.xlabel('stubby <-------------------> spiky (Z)')
				plt.ylabel('animate <-------------------> inanimate (Z)')
				plt.axvline(0, c='black')
				plt.axhline(0, c='black')
				plt.xlim(-3, 3)
				plt.ylim(-3, 3)
				plt.legend(frameon=False)#bbox_to_anchor=(0.9, 1.0))
				plt.tight_layout()
				plt.savefig(f'{resDir}/CatROIsPCspace.png')
				#plt.show()
				#plt.close()

				# nml regions
				preferences = {}

				resDir = os.path.join('analysis/results/fMRI/ROIsOnPCspace', topup, b0, HRFmodel, 'allSubjects')
				os.makedirs(resDir, exist_ok=True)
				fig, ax = plt.subplots(figsize=(4,4))

				for selectivity, colour in zip(['lat','med'], ['green', 'blue']):
					preferences[selectivity] = {}
					for r in range(2):
						region = f'nml_{selectivity}_{r+1}'
						shape = shapes[r]
						preferences[selectivity][region] = {}
						for component, contrast in zip(['animacy', 'aspectRatio'], ['animate>inanimate', 'spiky>stubby']):
							preferences[selectivity][region][component] = []
							for subject in experiment['scanInfo'].keys():
								theseData = dataAll[(dataAll['topup'] == topup) &
													(dataAll['b0'] == b0) &
													(dataAll['HRFmodel'] == HRFmodel) &
													(dataAll['nVoxels'] == voxelSize) &
													(dataAll['region'] == region) &
													(dataAll['condition'] == contrast) &
													(dataAll['subject'] == subject)]
								if len(theseData) > 0:
									if component == 'animacy':
										preferences[selectivity][region][component].append(-np.mean(theseData['PSC']))  # axis inverted
									elif component == 'aspectRatio':
										preferences[selectivity][region][component].append(np.mean(theseData['PSC']))
								else:
									preferences[selectivity][region][component].append(np.nan)

						ax.scatter(preferences[selectivity][region]['aspectRatio'], preferences[selectivity][region]['animacy'], s=30, color=colour, marker=shape)
				for selectivity, colour in zip(['lat', 'med'], ['green', 'blue']):
					for r in range(2):
						region = f'nml_{selectivity}_{r + 1}'
						shape = shapes[r]
						ax.scatter(np.nanmean(preferences[selectivity][region]['aspectRatio']), np.nanmean(preferences[selectivity][region]['animacy']), s=100,
								   edgecolors='black', marker=shape, color=colour, label=region)
				ax.set_title(f'NML ROIs in object space')
				plt.xlabel('stubby <-------------------> spiky (Z)')
				plt.ylabel('animate <-------------------> inanimate (Z)')
				plt.axvline(0, c='black')
				plt.axhline(0, c='black')
				plt.xlim(-3, 3)
				plt.ylim(-3, 3)
				plt.legend(frameon=False)#bbox_to_anchor=(0.9, .5))
				plt.tight_layout()
				plt.savefig(f'{resDir}/NML_PCspace.png')
				# plt.show()
				# plt.close()

				## category selective and NML regions

				resDir = os.path.join('analysis/results/fMRI/ROIsOnPCspace', topup, b0, HRFmodel, 'allSubjects')
				os.makedirs(resDir, exist_ok=True)
				fig, ax = plt.subplots(figsize=(5.4*0.85, 4*0.85))


				# nml
				preferences = {}

				for selectivity, colour in zip(['lat', 'med'], ['purple', 'orange']):
					preferences[selectivity] = {}
					for r in range(2):
						region = f'nml_{selectivity}_{r + 1}'
						regionLabel = f'NML {selectivity} {r + 1}'
						shape = shapes[r]
						preferences[selectivity][region] = {}
						for component, contrast in zip(['animacy', 'aspectRatio'], ['animate>inanimate', 'spiky>stubby']):
							preferences[selectivity][region][component] = []
							for subject in experiment['scanInfo'].keys():
								theseData = dataAll[(dataAll['topup'] == topup) &
													(dataAll['b0'] == b0) &
													(dataAll['HRFmodel'] == HRFmodel) &
													(dataAll['nVoxels'] == voxelSize) &
													(dataAll['region'] == region) &
													(dataAll['condition'] == contrast) &
													(dataAll['subject'] == subject)]
								if len(theseData) > 0:
									if component == 'animacy':
										preferences[selectivity][region][component].append(-np.mean(theseData['PSC']))  # axis inverted
									elif component == 'aspectRatio':
										preferences[selectivity][region][component].append(np.mean(theseData['PSC']))
								else:
									preferences[selectivity][region][component].append(np.nan)

						ax.scatter(preferences[selectivity][region]['aspectRatio'], preferences[selectivity][region]['animacy'], s=30, color=colour, marker=shape)

				preferencesNML = preferences

				# cat
				preferences = {}

				for selectivity, colour in zip(['body', 'face', 'house'], ['chartreuse', 'deepskyblue', 'red']):
					preferences[selectivity] = {}
					for r, region, in enumerate(selectivities[selectivity]):
						shape = shapes[r]
						preferences[selectivity][region] = {}
						for component, contrast in zip(['animacy', 'aspectRatio'], ['animate>inanimate', 'spiky>stubby']):
							preferences[selectivity][region][component] = []
							for subject in experiment['scanInfo'].keys():
								theseData = dataAll[(dataAll['topup'] == topup) &
													(dataAll['b0'] == b0) &
													(dataAll['HRFmodel'] == HRFmodel) &
													(dataAll['nVoxels'] == voxelSize) &
													(dataAll['region'] == region) &
													(dataAll['condition'] == contrast) &
													(dataAll['subject'] == subject)]
								if len(theseData) > 0:
									if component == 'animacy':
										preferences[selectivity][region][component].append(-np.mean(theseData['PSC']))  # axis inverted
									elif component == 'aspectRatio':
										preferences[selectivity][region][component].append(np.mean(theseData['PSC']))
								else:
									preferences[selectivity][region][component].append(np.nan)

						ax.scatter(preferences[selectivity][region]['aspectRatio'], preferences[selectivity][region]['animacy'], s=30, color=colour, marker=shape)

				preferencesCat = preferences
				# group means
				for selectivity, colour in zip(['lat', 'med'], ['purple','orange']):
					for r in range(2):
						region = f'nml_{selectivity}_{r + 1}'
						regionLabel = f'NML {selectivity} {r + 1}'
						shape = shapes[r]
						ax.scatter(np.nanmean(preferencesNML[selectivity][region]['aspectRatio']), np.nanmean(preferencesNML[selectivity][region]['animacy']), s=100,
								   edgecolors='black', marker=shape, color=colour, label=regionLabel)

				regionCounter=0
				for selectivity, colour in zip(['body', 'face', 'house'], ['chartreuse', 'deepskyblue', 'red']):
					for r, region, in enumerate(selectivities[selectivity]):
						shape = shapes[r]
						ax.scatter(np.nanmean(preferencesCat[selectivity][region]['aspectRatio']), np.nanmean(preferencesCat[selectivity][region]['animacy']), s=100,
								   edgecolors='black', marker=shape, color=colour, label=['EBA','FBA','OFA','FFA','pSTS','OPA','PPA','RSC'][regionCounter])
						regionCounter += 1
				ax.set_title(f'ROIs in object space')
				plt.xlabel('stubby <----------------> spiky (Z)')
				plt.ylabel('animate <-----------> inanimate (Z)')
				plt.axvline(0, c='black')
				plt.axhline(0, c='black')
				plt.xlim(-2, 2)
				plt.ylim(-2.7, 2.7)
				plt.legend(frameon=False,bbox_to_anchor=(1.0,1.01))  # bbox_to_anchor=(0.9, .5))
				plt.tight_layout()
				plt.savefig(f'{resDir}/NMLandCat_PCspace.png')
				# plt.show()
				# plt.close()


				# bar plot for category-selective regions
				resDir = os.path.join('analysis/results/fMRI/ROIresponses', topup, b0, HRFmodel, 'allSubjects')
				os.makedirs(resDir, exist_ok=True)
				regionOrder = ['eba','fba','ofa','ffa','psts','opa','ppa','rsc']
				labels = ['EBA','FBA','OFA','FFA','pSTS','OPA','PPA','RSC']
				for scan in experiment['design'].keys():
					if scan == 'animacyAspectRatioLoc_v2':
						conds = ['animate_spiky','animate_stubby','inanimate_spiky','inanimate_stubby']
						condsLabels = ['spiky-animate','stubby-animate','spiky-inanimate','stubby-inanimate']
						colours = ['chartreuse', 'deepskyblue', 'orange','fuchsia']
					elif scan == 'BFHOloc_v2':
						conds = ['body', 'face','house','object']
						condsLabels=conds
						colours = ['chartreuse', 'deepskyblue', 'red', 'gold']

					df = dataAll[(dataAll['topup'] == topup) &
								 (dataAll['b0'] == b0) &
								 (dataAll['HRFmodel'] == HRFmodel) &
								 (dataAll['nVoxels'] == voxelSize) &
								 (dataAll['region'].isin(regionOrder)) &
								 (dataAll['condition'].isin(conds))]
					df = df.groupby(['subject', 'region', 'condition']).agg(['mean']).reset_index().iloc[:,[0,1,2,5]]
					df = df.groupby(['region', 'condition']).agg(['mean','sem']).reset_index()
					df.columns = df.columns[:2].droplevel([1,2]).append(df.columns[2:4].droplevel([0,1]))
					dfMeans = df.pivot(index='region', columns='condition', values='mean').reindex(index = regionOrder)
					dfSEMs = df.pivot(index='region', columns='condition', values='sem').reindex(index = regionOrder).values
					dfPlot = dfMeans.plot(kind='bar', color=colours, yerr=dfSEMs.transpose(), ylabel='signal change (%)', rot=0, figsize = (11*0.7,4*0.7))
					fig = dfPlot.get_figure()
					plt.ylabel('signal change (%)')
					plt.xlabel(None)
					plt.legend(frameon=False, bbox_to_anchor=(1.01,1), labels = condsLabels).set_title('')
					plt.xticks(np.arange(len(regionOrder)), labels=labels)
					plt.ylim(-1,5)
					plt.tight_layout()
					#plt.gca().invert_yaxis()
					fig.savefig(os.path.join(resDir, f'{scan}_barPlot.png'))
					#plt.show(block=False)
					#plt.close()

				# bar plot for nml regions
				regionOrder = selectivities['spiky']
				for scan in experiment['design'].keys():
					if scan == 'animacyAspectRatioLoc_v2':
						conds = ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']
						condsLabels = ['spiky-animate', 'stubby-animate', 'spiky-inanimate', 'stubby-inanimate']
						colours = ['chartreuse', 'deepskyblue', 'orange', 'fuchsia']
					elif scan == 'BFHOloc_v2':
						conds = ['body', 'face', 'house', 'object']
						condsLabels = conds
						colours = ['chartreuse', 'deepskyblue', 'red', 'gold']

					df = dataAll[(dataAll['topup'] == topup) &
								 (dataAll['b0'] == b0) &
								 (dataAll['HRFmodel'] == HRFmodel) &
								 (dataAll['nVoxels'] == voxelSize) &
								 (dataAll['region'].isin(regionOrder)) &
								 (dataAll['condition'].isin(conds))]
					df = df.groupby(['subject', 'region', 'condition']).agg(['mean']).reset_index().iloc[:, [0, 1, 2, 5]]
					df = df.groupby(['region', 'condition']).agg(['mean', 'sem']).reset_index()
					df.columns = df.columns[:2].droplevel([1, 2]).append(df.columns[2:4].droplevel([0, 1]))
					dfMeans = df.pivot(index='region', columns='condition', values='mean').reindex(index=regionOrder)
					dfSEMs = df.pivot(index='region', columns='condition', values='sem').reindex(index=regionOrder).values
					dfPlot = dfMeans.plot(kind='bar', color=colours, yerr=dfSEMs.transpose(), ylabel='signal change (%)', rot=0, figsize=(7, 4))
					fig = dfPlot.get_figure()
					plt.ylabel('signal change (%)')
					plt.legend(frameon=False, labels = condsLabels, bbox_to_anchor=(1.0, 1.01)).set_title('')
					# plt.gca().invert_yaxis()
					fig.savefig(os.path.join(resDir, f'{scan}_barPlot_NML.png'))
				# plt.show(block=False)
				# plt.close()

				# bar plot for nml regions without category selective voxels
				regionOrder = selectivities['spiky'].copy()
				for r, regionO in enumerate(regionOrder):
					regionOrder[r] = regionO + '_noCatVoxels'
				for scan in experiment['design'].keys():
					if scan == 'animacyAspectRatioLoc_v2':
						conds = ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']
						colours = ['chartreuse', 'deepskyblue', 'orange', 'fuchsia']
					elif scan == 'BFHOloc_v2':
						conds = ['body', 'face', 'house', 'object']
						colours = ['chartreuse', 'deepskyblue', 'red', 'gold']

					df = dataAll[(dataAll['topup'] == topup) &
								 (dataAll['b0'] == b0) &
								 (dataAll['HRFmodel'] == HRFmodel) &
								 (dataAll['nVoxels'] == voxelSize) &
								 (dataAll['region'].isin(regionOrder)) &
								 (dataAll['condition'].isin(conds))]
					# add mean of all regions
					dfAllNML = df.groupby(['subject', 'condition']).agg(['mean']).reset_index().iloc[:, [0, 1, 4]]
					dfAllNML['region'] = 'NML'
					df = df.groupby(['subject', 'region', 'condition']).agg(['mean']).reset_index().iloc[:, [0, 1, 2, 5]]
					df = pd.concat([dfAllNML, df])
					df = df.groupby(['region', 'condition']).agg(['mean', 'sem']).reset_index()
					df.columns = df.columns[:2].droplevel([1, 2]).append(df.columns[2:4].droplevel([0, 1]))
					dfMeans = df.pivot(index='region', columns='condition', values='mean').reindex(index=['NML'] + regionOrder)
					dfSEMs = df.pivot(index='region', columns='condition', values='sem').reindex(index=['NML'] + regionOrder).values
					dfPlot = dfMeans.plot(kind='bar', color=colours, yerr=dfSEMs.transpose(), ylabel='signal change (%)', rot=0, figsize=(7*0.7, 4*0.7))
					fig = dfPlot.get_figure()
					plt.ylabel('signal change (%)')
					plt.xlabel(None)
					plt.ylim(0, 5)
					plt.xticks(ticks = np.arange(5),labels = ['NML\nall regions', 'NML\nlat 1','NML\nlat 2','NML\nmed 1','NML\nmed 2'])

					plt.legend(frameon=False, bbox_to_anchor=(1.01, 1)).set_title('')

					plt.tight_layout()
					# plt.gca().invert_yaxis()
					fig.savefig(os.path.join(resDir, f'{scan}_barPlot_NMLnoCat.png'))
			# plt.show(block=False)
			# plt.close()
