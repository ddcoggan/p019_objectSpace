''''

'''

import os
import sys
import glob
import datetime
import shutil
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt


def calculateOverlapROI(pathA, pathB):

	dataA = nib.load(pathA).get_fdata().flatten()
	dataB = nib.load(pathB).get_fdata().flatten()
	nVoxA = np.count_nonzero(dataA)
	nVoxB = np.count_nonzero(dataB)
	overlap = sum(dataA * dataB)
	return(nVoxA, nVoxB, overlap)


# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

sizes = [16, 64, 256]
thr = 3.1
stdMaskDir = '/mnt/HDD12TB/masks'
VSpath = '/mnt/HDD12TB/masks/HOCSA/ventral_stream.nii.gz'
regions = {'BFHOloc_v2': {'body': ['eba', 'fba'],
		   				  'face': ['ofa', 'ffa', 'psts'],
		   				  'house': ['opa', 'ppa', 'rsc']},
		   'animacyAspectRatioLoc_v2': {'NML': ['nml_lat_1', 'nml_lat_2', 'nml_med_1', 'nml_med_2']}}

topup = 'topUp'
b0 = 'noB0'
HRFmodel = 'doubleGamma'

# KNOWN BFHO REGIONS VERSUS AAR CLUSTERS
print('Analysing known BFHO selective ROIs v whole-brain AAR clusters...')
reportFile = 'data/fMRI/group/overlapROIknownRegions.csv'
mi = open(reportFile, 'w+')
mi.write('subject,session,selectivity,region,hemi,nVoxTarget,condition,voxType,nVox\n')

for subject in experiment['scanInfo'].keys():
	for s, session in enumerate(experiment['scanInfo'][subject].keys()):

		# just use first run directory
		sessDir = os.path.join('data/fMRI/individual', subject, session)

		for category in regions['BFHOloc_v2']:
			for region in regions['BFHOloc_v2'][category]:
				for hemi in ['lh','rh']:

					# do this for final ROIs
					for targetSize in [16, 64, 256]:
						regionPath = os.path.join(sessDir, f'masks/native/{region}_{hemi}_{targetSize:05}_voxels.nii.gz')
						if os.path.isfile(regionPath):

							for c, cond in enumerate(['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']):

								condPath = os.path.join(sessDir, f'masks/native/{cond}Selective.nii.gz')
								nVoxRegion, nVoxCond, overlap = calculateOverlapROI(regionPath, condPath)
								mi.write(f'{subject},{session},{category},{region},{hemi},{targetSize},{cond},regionUnique,{nVoxRegion-overlap}\n')
								mi.write(f'{subject},{session},{category},{region},{hemi},{targetSize},{cond},overlap,{overlap}\n')
								mi.write(f'{subject},{session},{category},{region},{hemi},{targetSize},{cond},conditionUnique,{nVoxCond-overlap}\n')

					# and for entire ROI activation
					maskedActPath = os.path.join(sessDir, f'masks/native/floodfill/masked_activations/{region}_{hemi}.nii.gz')
					if os.path.isfile(maskedActPath):
						regionPath = f'{maskedActPath[:-7]}_thr_bin.nii.gz'
						if not os.path.isfile(regionPath):
							os.system(f'fslmaths {maskedActPath} -thr {thr} -bin {regionPath}')

						for c, cond in enumerate(['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']):
							condPath = os.path.join(sessDir, f'masks/native/{cond}Selective.nii.gz')
							nVoxRegion, nVoxCond, overlap = calculateOverlapROI(regionPath, condPath)
							mi.write(f'{subject},{session},{category},{region},{hemi},allVoxels,{cond},regionUnique,{nVoxRegion - overlap}\n')
							mi.write(f'{subject},{session},{category},{region},{hemi},allVoxels,{cond},overlap,{overlap}\n')
							mi.write(f'{subject},{session},{category},{region},{hemi},allVoxels,{cond},conditionUnique,{nVoxCond - overlap}\n')
mi.close()

# NML REGIONS VERSUS BFHO CLUSTERS
print('Analysing NML ROIs v whole-brain BFHO clusters...')
reportFile = 'data/fMRI/group/overlapROINMLRegions.csv'
mi = open(reportFile, 'w+')
mi.write('subject,session,selectivity,region,hemi,nVoxTarget,condition,nVox\n')

for subject in experiment['scanInfo'].keys():
	for s, session in enumerate(experiment['scanInfo'][subject].keys()):

		# just use first run directory
		sessDir = os.path.join('data/fMRI/individual', subject, session)


		for category in regions['animacyAspectRatioLoc_v2']:
			for region in regions['animacyAspectRatioLoc_v2'][category]:
				for hemi in ['lh','rh']:

					# do this for final ROIs
					for targetSize in [16, 64, 256]:
						regionPath = os.path.join(sessDir, f'masks/native/{region}_{hemi}_{targetSize:05}_voxels.nii.gz')
						if os.path.isfile(regionPath):

							totalOverlap = 0
							for c, cond in enumerate(['body','face','house']):

								condPath = os.path.join(sessDir, f'masks/native/{cond}Selective.nii.gz')
								nVoxRegion, nVoxCond, overlap = calculateOverlapROI(regionPath, condPath)
								mi.write(f'{subject},{session},{category},{region},{hemi},{targetSize},{cond},{overlap}\n')
								totalOverlap += overlap
							mi.write(f'{subject},{session},{category},{region},{hemi},{targetSize},unique,{nVoxRegion-totalOverlap}\n')

					# and for entire ROI activation
					maskedActPath = os.path.join(sessDir, f'masks/native/floodfill/masked_activations/{region}_{hemi}.nii.gz')
					if os.path.isfile(maskedActPath):
						regionPath = f'{maskedActPath[:-7]}_thr_bin.nii.gz'
						if not os.path.isfile(regionPath):
							os.system(f'fslmaths {maskedActPath} -thr {thr} -bin {regionPath}')

						totalOverlap = 0
						for c, cond in enumerate(['body','face','house']):

							condPath = os.path.join(sessDir, f'masks/native/{cond}Selective.nii.gz')
							nVoxRegion, nVoxCond, overlap = calculateOverlapROI(regionPath, condPath)
							mi.write(f'{subject},{session},{category},{region},{hemi},allVoxels,{cond},{overlap}\n')
							totalOverlap += overlap
						mi.write(f'{subject},{session},{category},{region},{hemi},allVoxels,unique,{nVoxRegion - totalOverlap}\n')

# BFHO CLUSTERS V AAR CLUSTERS
print('Analysing known BFHO cluster v AAR clusters in ventral stream...')
reportFile = 'data/fMRI/group/overlapROIallClusters.csv'
mi = open(reportFile, 'w+')
mi.write('subject,session,BFHOcond,AARcond,voxType,nVox\n')

for subject in experiment['scanInfo'].keys():
	for s, session in enumerate(experiment['scanInfo'][subject].keys()):

		# just use first run directory
		sessDir = os.path.join('data/fMRI/individual', subject, session)


		for category in ['body', 'face', 'house', 'object']:
			regionPath = os.path.join(sessDir,f'masks/native/{category}SelectiveVS.nii.gz')

			for c, cond in enumerate(['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']):
				condPath = os.path.join(sessDir, f'masks/native/{cond}SelectiveVS.nii.gz')
				nVoxRegion, nVoxCond, overlap = calculateOverlapROI(regionPath,
																	condPath)
				mi.write(f'{subject},{session},{category},{cond},BFHOunique,{nVoxRegion-overlap}\n')
				mi.write(f'{subject},{session},{category},{cond},overlap,{overlap}\n')
				mi.write(f'{subject},{session},{category},{cond},AARunique,{nVoxCond - overlap}\n')
mi.close()



