#!/usr/bin/python
'''
collates featquery reports. is done separately from featquery analyses to allow all results to be collected at once, not just from analyses actually run at the time.
this should be fixable in the previous script but doesnt make a difference really
'''

import os
import glob
import shutil
import numpy as np

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

resultsFile = 'data/fMRI/group/featquery.csv'

results = open(resultsFile, 'w')
results.write('topup,b0,HRFmodel,subject,session,scan,run,roi,condition,PSC,lower,upper\n')

for topup in ['topUp']:#, 'noTopUp']:
	for b0 in ['noB0']: #, 'b0']:
		for HRFmodel in ['doubleGamma']: #, 'singleGamma']:
			for subject in experiment['scanInfo'].keys():
				for session in experiment['scanInfo'][subject].keys():
					sessDir = os.path.join('data/fMRI/individual', subject, session)
					masks = sorted(glob.glob(os.path.join(sessDir, 'masks/native/*64*')))  # go by masks rather than featquery dirs to make sure all featqueries have run

					for scan in experiment['design'].keys():

						if scan == 'BFHOloc_v2':
							condNames = ['body', 'face', 'house', 'object', 'face>body', 'face>house', 'face>object',
										 'allConds', 'body>others', 'face>others', 'house>others', 'object>others']
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
								reportFile = os.path.join(runDir, 'featquery', maskName, 'report.txt')
								if os.path.isfile(reportFile):
									report = open(reportFile, 'r')
									report = report.readlines()
									for c, cond in enumerate(condNames):
										info = report[c].split(' ')
										lower, PSC, upper = np.array(info)[[4,5,7]]
										results.write(f'{topup},{b0},{HRFmodel},{subject},{session},{scan},{run+1},{maskName},{cond},{PSC},{lower},{upper}\n')
results.close()

