#!/usr/bin/python
'''
takes first level design files (each run) made prior to this script, edits and submits to feat
'''

import os
import datetime
import time
import shutil

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

# first level analysis
for topup, topupString in zip(['topUp'],['_topUp']): #zip(['topUp', 'noTopUp'],['_topUp','']):
    for b0, b0String in zip(['noB0'],['_b0']): #zip(['b0', 'noB0'],['_b0', '']):
        for HRFmodel in ['doubleGamma']: #, 'singleGamma']:
            for subject in experiment['scanInfo'].keys():
                for session in experiment['scanInfo'][subject].keys():
                    for scan in experiment['design'].keys():
                        if scan != 'restingState':
                            for r, run in enumerate(experiment['scanInfo'][subject][session]['funcScans'][scan]):

                                print(f'\n{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | Run: {r+1} | Top Up: {topup} |  B0: {b0} | HRF model: {HRFmodel}')

                                outDir = os.path.join('data/fMRI/individual', subject, session, 'functional', scan, f'run{r+1:02}/outputs', topup, b0, HRFmodel, 'firstLevel.feat')

                                if not os.path.isdir(outDir):
                                    print('Analysis not found or was just deleted, analysing...')

                                    # replace relevant settings in design file. The original should be set up for the first run of the first subject,
                                    # and the files should be set up so that all that needs changing is the subject ID and the run number
                                    designFile = os.path.join('data/fMRI/designs/firstLevel', scan, HRFmodel, 'design.fsf')
                                    with open(designFile, 'r') as file:
                                        fileData = file.read()

                                    # replace topup and b0 types for input field and output field
                                    fileData = fileData.replace('rawData', f'rawData{topupString}{b0String}')
                                    fileData = fileData.replace('noTopUp', topup)
                                    fileData = fileData.replace('noB0', b0)

                                    # replace subject ID
                                    fileData = fileData.replace('F117', subject)

                                    # replace session
                                    fileData = fileData.replace('210223', session)

                                    # replace run number
                                    fileData = fileData.replace('run01', f'run{r+1:02}')

                                    # check out dir in design file is correct
                                    lines = fileData.splitlines()
                                    for line in lines:
                                        if line.startswith(f'set fmri(outputdir)'):
                                            actualOutDir = line.split()[2][1:-1]
                                            if not actualOutDir == os.path.join(os.getcwd(), outDir):
                                                print(os.path.join(os.getcwd(), outDir))
                                                print(actualOutDir)
                                                raise Exception(f'Destination directory does not match, check feat analysis at {outDir}')

                                    # write the file out again
                                    designFileTemp = f'{designFile[:-4]}_temp.fsf'
                                    with open(designFileTemp, 'w') as file:
                                        file.write(fileData)
                                    file.close()

                                    # run analysis
                                    os.system(f'feat {designFileTemp}')

                                    # check that analysis has completed successfully
                                    if not os.path.isfile(os.path.join(outDir, 'stats/cope1.nii.gz')):
                                        raise Exception(f'No cope files produced, check feat analysis at {outDir}')


                                elif time.time() - os.stat(outDir).st_ctime > 1800:  # if created over 2 hours ago

                                        # check that analysis has completed successfully

                                        if not os.path.isfile(os.path.join(outDir, 'stats/cope1.nii.gz')):

                                            raise Exception(f'No cope files found in previously completed feat dir, check feat analysis at {outDir}')

                                        else:

                                            print('Analysis with cope files found, skipping...')
                                else:

                                    print('analysis started < 30 mins ago, not checked')

print('Done')
