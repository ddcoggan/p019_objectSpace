'''
makes flood-filled masks and other masks
script runs from scratch each time, overwriting all previous files. It needs to do this in order to collate statistics of all masks in report files.
'''

import os
import sys
import glob
import datetime
import shutil
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('/mnt/HDD12TB/masterScripts/fMRI')
from splitByHemi import splitByHemi
from makeFloodfillMasks import makeFloodfillMasks

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

overwrite = False
sizes = [16, 64, 256]
thr = 3.1
stdMaskDir = '/mnt/HDD12TB/masks'
VSpath = '/mnt/HDD12TB/masks/HOCSA/ventral_stream.nii.gz'
locInfo = {'BFHOloc_v2': {'body': {'regions': ['eba', 'fba'],
								'locCope': 9},
		   			   'face': {'regions': ['ofa', 'ffa', 'psts'],
								'locCope': 10},
		   			   'house': {'regions': ['opa', 'ppa', 'rsc'],
								 'locCope': 11}},
		   'animacyAspectRatioLoc_v2': {'NML': {'regions': ['nml_lat_1', 'nml_lat_2', 'nml_lat_3', 'nml_med_1', 'nml_med_2', 'nml_med_3'],
												  'locCope': 14}}}


for topup in ['topUp']:#, 'noTopUp']:
	for b0 in ['noB0']:#, 'b0']:
		for HRFmodel in ['doubleGamma']: #, 'singleGamma']:
			for subject in experiment['scanInfo']:
				for s, session in enumerate(experiment['scanInfo'][subject].keys()):

					# key directories
					sessDir = os.path.join('data/fMRI/individual', subject, session)
					stdMaskDir = f'{sessDir}/masks/standard'
					os.makedirs(stdMaskDir, exist_ok=True)
					estimateDir = f'{stdMaskDir}/floodfill/estimates'
					natMaskDir = f'{sessDir}/masks/native'
					os.makedirs(natMaskDir, exist_ok=True)
					floodfillDir = os.path.join(natMaskDir, 'floodfill')
					os.makedirs(floodfillDir, exist_ok=True)
					refFunc = os.path.join(sessDir, 'reg/refFunc.nii.gz')

					# record final mask stats
					reportFile = os.path.join(floodfillDir, 'maskInfo.csv')
					mi = open(reportFile, 'w+')
					mi.write('region,hemi,peakX,peakY,peakZ,peakX(MNI),peakY(MNI),peakZ(MNI),Nvox(target),Nvox(actual),minZ\n')  # write header (actual data is added by floodfillMasks.py)
					mi.close()

					# make non-floodfill masks based on significant activation in ventral stream
					for scan in experiment['design'].keys():
						if scan == 'BFHOloc_v2':
							condNames = ['body', 'face', 'house', 'object']
						elif scan == 'animacyAspectRatioLoc_v2':
							condNames = ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby', 'NML']

						for c, cond in enumerate(condNames):

							if cond == 'NML':
								buffer = 10
							else:
								buffer = 9

							# get cluster correction mask an binarize to remove cluster IDs
							clusterMask = os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'allRuns', topup, b0, HRFmodel, f'{thr}/secondLevel.gfeat/cope{c + buffer}.feat/cluster_mask_zstat1.nii.gz')
							clusterMaskBin = f'{os.path.dirname(clusterMask)}/cluster_mask_zstat1_bin.nii.gz'
							if not os.path.isfile(clusterMaskBin):
								os.system(f'fslmaths {clusterMask} -bin {clusterMaskBin}')

							# get selective masks across whole brain
							outPath = f'{stdMaskDir}/{cond}Selective.nii.gz'
							if not os.path.exists(outPath):
								shutil.copy(clusterMaskBin, outPath)

							# get selective masks within ventral stream
							outPath = f'{stdMaskDir}/{cond}SelectiveVS.nii.gz'
							if not os.path.exists(outPath):
								os.system(f'fslmaths {clusterMaskBin} -mul {VSpath} {outPath}')


						# native space masks including flood fill
						# get sample run directory
						runDir = os.path.join('data/fMRI/individual', subject, session, 'functional', scan,
											  f'run01/outputs', topup, b0, HRFmodel, 'firstLevel.feat')
						transMat = os.path.join(runDir, 'reg/standard2example_func.mat')
						
						# for ROI plots
						EFrange = os.popen(f'fslstats {refFunc} -R')
						EFmax = float(EFrange.read().split()[1])

						# convert all non-floodfill masks to native space
						for cond in condNames:

							# ventral stream
							roiPathStd = f'{stdMaskDir}/{cond}SelectiveVS.nii.gz'
							roiPathNat = os.path.join(natMaskDir, f'{cond}SelectiveVS.nii.gz')
							if not os.path.isfile(roiPathNat):
								print(f'converting {cond}SelectiveVS to native space')
								os.system(f'flirt -in {roiPathStd} -ref {refFunc} -out {roiPathNat} -init {transMat} -applyxfm -interp nearestneighbour')

							# whole brain
							roiPathStd = f'{stdMaskDir}/{cond}Selective.nii.gz'
							roiPathNat = os.path.join(natMaskDir, f'{cond}Selective.nii.gz')
							print(f'converting {cond}Selective to native space')
							os.system(f'flirt -in {roiPathStd} -ref {refFunc} -out {roiPathNat} -init {transMat} -applyxfm -interp nearestneighbour')

						# floodfill masks
						for condName in locInfo[scan]:

							# convert activation map to native func
							actPathStd = os.path.join(sessDir, f'functional/{scan}/allRuns/{topup}/{b0}/{HRFmodel}/{thr}/secondLevel.gfeat/cope{locInfo[scan][condName]["locCope"]}.feat/stats/zstat1.nii.gz')
							actPathNat = os.path.join(floodfillDir, f'actMap_{condName}.nii.gz')
							print(f'converting actMap_{condName} to native space')
							os.system(f'flirt -in {actPathStd} -ref {refFunc} -out {actPathNat} -init {transMat} -applyxfm -interp trilinear')

							# convert clusterMask to native func
							clusterMaskStd = os.path.join(sessDir, f'functional/{scan}/allRuns/{topup}/{b0}/{HRFmodel}/{thr}/secondLevel.gfeat/cope{locInfo[scan][condName]["locCope"]}.feat/cluster_mask_zstat1_bin.nii.gz')
							if not os.path.exists(clusterMaskStd):
								os.system(f'fslmaths {os.path.dirname(clusterMaskStd)}/cluster_mask_zstat1.nii.gz -bin {clusterMaskStd}')
							clusterMaskNat =  os.path.join(floodfillDir, f'clusterMask_{condName}.nii.gz')
							print(f'converting clusterMask_{condName} to native space')
							os.system(f'flirt -in {clusterMaskStd} -ref {refFunc} -out {clusterMaskNat} -init {transMat} -applyxfm -interp nearestneighbour')

							# combine
							actPathNatCC = os.path.join(floodfillDir, f'actMap_{condName}_CC.nii.gz')
							os.system(f'fslmaths {actPathNat} -mul {clusterMaskNat} {actPathNatCC}')

							for region in locInfo[scan][condName]['regions']:

								for hemi in ['lh','rh']:

									estimateStd = os.path.join(estimateDir, f'{region}_{hemi}.nii.gz')

									if os.path.isfile(estimateStd):

										# get peak MNI coords for each mask
										coordsStd = os.popen(f'fslstats {actPathStd} -k {estimateStd} -x')
										coordsStd = coordsStd.read()[:-2].split(' ')

										# convert estimate to native space
										estimateNat = os.path.join(floodfillDir, f'{region}_{hemi}.nii.gz')
										print(f'converting {region} {hemi} estimate to native space')
										os.system(f'flirt -in {estimateStd} -ref {refFunc} -out {estimateNat} -init {transMat} -applyxfm -interp nearestneighbour')

										# mask the native activation by the native ROI estimate
										maskedActivationDir = os.path.join(floodfillDir, 'masked_activations')
										os.makedirs(maskedActivationDir, exist_ok=True)
										maskedActivation = os.path.join(maskedActivationDir, f'{region}_{hemi}.nii.gz')
										print(f'converting {region} {hemi} masked activation to native space')
										os.system(f'fslmaths {actPathNatCC} -mul {estimateNat} {maskedActivation}')

										# make distribution plot of t values in masked activation
										actData = nib.load(maskedActivation).get_fdata().flatten()
										actData.sort()

										nVoxels = np.count_nonzero(np.maximum(0,actData)) # reLu then sum non zero
										y = actData[:-(nVoxels+1):-1]
										x = np.arange(1,nVoxels+1)

										figPath = os.path.join(maskedActivationDir, f'{region}_{hemi}_zvals.png')
										plt.plot(x, y)
										plt.xlabel('voxels')
										plt.ylabel('t value')
										plt.title(f'subject: {subject}, session: {s+1}. region: {region}, hemisphere: {hemi}')
										plt.grid(True)
										plt.savefig(figPath)
										plt.show(block=False)
										plt.close()

										# get location of peak activation
										coords = os.popen(f'fslstats {maskedActivation} -x')
										coords = [int(x) for x in coords.read()[:-2].split(' ')]

										# run flood-fill algorithm
										for size in sizes:

											maskOutFile = os.path.join(natMaskDir, f'{region}_{hemi}_{size:05}_voxels.nii.gz')
											print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Top Up: {topup} |  B0: {b0} | '
												  f'HRF model: {HRFmodel} | Subject: {subject} | Session: {session} | Scan: {scan} | '
												  f'Condition: {condName} | Region: {region} | Hemi: {hemi} | Size: {size}')

											if not os.path.isfile(maskOutFile) or overwrite:
												sizeAct, minZ = makeFloodfillMasks(maskedActivation,size,coords,thr,maskOutFile)
											else:
												sizeAct = os.popen(f'fslstats {maskOutFile} -V').read().split(' ')[0]
												minZ = os.popen(f'fslstats {maskedActivation} -k {maskOutFile} -R').read().split(' ')[0]

											# record mask info
											mi = open(reportFile, 'a')
											mi.write(f'{region},{hemi},{coords[0]},{coords[1]},{coords[2]},{coordsStd[0]},{coordsStd[1]},{coordsStd[2]},{size},{sizeAct},{minZ}\n')
											mi.close()

											# make brain plots showing ROI
											plotDir = os.path.join(floodfillDir, f'plots')
											os.makedirs(plotDir, exist_ok=True)
											plotFile = os.path.join(plotDir, f'{region}_{hemi}_{size:05}_voxels.png')
											maxAct = os.popen(f'fslstats {actPathNat} -R')
											maxAct = float(maxAct.read().split()[1])


											fsleyesCommand = f'fsleyes render --outfile {plotFile} --size 3200 600 --scene ortho --autoDisplay -vl {coords[0]} {coords[1]} {coords[2]} ' \
														 f'{refFunc} -dr 0 {EFmax} {actPathNat} -dr {thr} {maxAct} -cm red-yellow {maskedActivation} ' \
														 f'-dr {thr} {maxAct} -cm blue-lightblue {maskOutFile} -dr 0 1 -cm green'
											os.system(fsleyesCommand)

					# make NML masks of voxels not significant for classic categories
					NMLmasks = glob.glob(f'{natMaskDir}/nml*')
					for NMLmask in NMLmasks:
						if not NMLmask.endswith('noCatVoxels.nii.gz'):
							outFile = NMLmask[:-7]+'_noCatVoxels.nii.gz'
							os.system(f'fslmaths {NMLmask} -sub {natMaskDir}/bodySelective.nii.gz -sub {natMaskDir}/faceSelective.nii.gz -sub {natMaskDir}/houseSelective.nii.gz -thr 1 {outFile}')

print('FINISHED')


