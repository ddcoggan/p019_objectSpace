#!/usr/bin/python
'''
runs featquery analysis
'''

import os
import glob
import shutil
import datetime
import time
#time.sleep(3600)

overwrite = False

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

for topup in ['topUp']: #, 'noTopUp']:
	for b0 in ['noB0']: #, 'b0']:
		for HRFmodel in ['doubleGamma']: #, 'singleGamma']:
			for subject in experiment['scanInfo']:
				for session in experiment['scanInfo'][subject]:
					sessDir = os.path.join('data/fMRI/individual', subject, session)
					masks = sorted(glob.glob(f'{sessDir}/masks/native/*.nii.gz'))

					for scan in experiment['design']:
						if scan == 'BFHOloc_v2':
							condNames = ['body', 'face', 'house', 'object',
                                     'face>body', 'face>house', 'face>object', 'allConds',
                                     'body>others', 'face>others', 'house>others', 'object>others']
						elif scan == 'animacyAspectRatioLoc_v2':
							condNames = ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby',
                                                   'inanimate_spiky>inanimate_stubby', 'animate>inanimate', 'spiky>stubby', 'allConds',
                                                   'animate_spiky>others', 'animate_stubby>others', 'inanimate_spiky>others', 'inanimate_stubby>others',
                                                   'animate_spiky>animate_stubby', 'NML',
                                                   'animate_stubby>inanimate_spiky','inanimate_stubby>animate_spiky','inanimate_spiky>animate_stubby', 'animate_spiky>inanimate_stubby']

						for run in range(len(experiment['scanInfo'][subject][session]['funcScans'][scan])):
							runDir = os.path.join('data/fMRI/individual', subject, session, 'functional', scan,
												  f'run{run + 1:02}/outputs', topup, b0, HRFmodel, 'firstLevel.feat')


							for mask in masks:
								maskName = os.path.basename(mask)[:-7]
								print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | topup: {topup} | B0: {b0} | HRFmodel: {HRFmodel} | Subject: {subject} | Session: {session} | Scan: {scan} | Run: {run+1} | Region: {maskName} ')
								if not os.path.exists(os.path.join(runDir, 'featquery')):
									os.makedirs(os.path.join(runDir, 'featquery'), exist_ok=True)

								outDir = os.path.join(runDir, 'featquery', maskName)
								if os.path.isdir(outDir) and overwrite:
									print('Deleting previous analysis directory...')
									shutil.rmtree(outDir)

								nCopes = len(glob.glob(os.path.join(runDir, 'stats/cope*.nii.gz')))
								fqCommand = f'featquery 1 {os.path.join(os.getcwd(),runDir)} {nCopes}'
								for cope in range(nCopes):
									fqCommand += f' stats/cope{cope+1}'
								fqCommand += f' featquery/{maskName} -p {os.path.join(os.getcwd(), sessDir, "masks/native", maskName)}'

								if not os.path.isdir(outDir) or overwrite:
									os.system(fqCommand)

print('Done.')



