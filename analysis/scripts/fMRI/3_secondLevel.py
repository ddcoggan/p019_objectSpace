
'''
takes second level design files (combining analyses across runs) made prior to this script, edits and submits to feat
'''

import os
import glob
import pickle
import datetime
import shutil
import time
#time.sleep(3600)

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

for subject in experiment['scanInfo'].keys():
    for session in experiment['scanInfo'][subject].keys():
        for scan in experiment['design'].keys():
            for topup in ['topUp']: #['topUp', 'noTopUp']
                for b0 in ['noB0']: # ['b0', 'noB0']:
                    for HRFmodel in ['doubleGamma']: # ['doubleGamma', 'singleGamma']:
                        for thr in ['2.3','2.7','3.1','3.5']:

                            print(f'\n{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | topup: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                            # get list of runs for first subject (only used to identify string patterns in design file to replace later)
                            # Note: second session is used as first session for M012 contained one run of each scan
                            runDirsTemplate = sorted(glob.glob(os.path.join('data/fMRI/individual/F117/210223/functional', scan, 'run??/outputs/noTopUp/noB0/doubleGamma/firstLevel.feat')))
                            nRunsTemplate = len(runDirsTemplate)

                            # compile list of runs for this subject
                            runDirs = sorted(glob.glob(os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'run??/outputs', topup, b0, HRFmodel, 'firstLevel.feat')))
                            nRuns = len(runDirs)

                            if nRuns > 1:

                                outDir = os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'allRuns', topup, b0, HRFmodel, thr, 'secondLevel.gfeat')
                                os.makedirs(os.path.dirname(outDir), exist_ok=True)

                                if not os.path.isdir(outDir):

                                    print('Analysis not found, analysing...')
                                    designFile = os.path.join('data/fMRI/designs/secondLevel', scan, 'design.fsf')

                                    with open(designFile, 'r') as file:
                                        fileData = file.read()

                                    # replace subject ID
                                    fileData = fileData.replace('F117', subject)

                                    # replace session
                                    fileData = fileData.replace('210223', session)

                                    # replace topup and b0 types for input field and output field
                                    fileData = fileData.replace('topUp', topup)
                                    fileData = fileData.replace('noB0', b0)

                                    # replace HRFmodel
                                    fileData = fileData.replace('doubleGamma', HRFmodel)

                                    # replace threshold
                                    fileData = fileData.replace('set fmri(z_thresh) 3.1', f'set fMRI(z_thresh) {thr}')
                                    fileData = fileData.replace('/3.1/', f'/{thr}/')

                                    # replace number of runs
                                    fileData = fileData.replace(f'set fmri(npts) {nRunsTemplate}', f'set fmri(npts) {nRuns}')
                                    fileData = fileData.replace(f'set fmri(multiple) {nRunsTemplate}', f'set fmri(multiple) {nRuns}')

                                    # add/remove input analysis directories if different to template
                                    if len(runDirs) > len(runDirsTemplate):
                                        for r, runDir in enumerate(runDirs):
                                            fileData += f'\nset fmri(evg{r + 1}.1) 1.0'
                                            fileData += f'\nset fmri(groupmem.{r + 1}) 1'
                                            fileData += f'\nset feat_files({r + 1}) "{os.path.join(os.getcwd(), runDir)}"'

                                    elif nRuns < nRunsTemplate:
                                        for r in range(nRuns, nRunsTemplate):
                                            fileData = fileData.replace(f'set feat_files({r+1})', '#')
                                            fileData = fileData.replace(f'set fmri(evg{r+1}.', '#')
                                            fileData = fileData.replace(f'set fmri(groupmem.{r + 1})', '#')

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
                                    if not os.path.isfile(os.path.join(outDir, 'cope1.feat/stats/cope1.nii.gz')):
                                        raise Exception(f'No cope files produced, check feat analysis at {outDir}')

                            else:
                                print('Only a single run found for this session, skipping...')

print('Done')
