'''
takes third level design files (combining analyses across subjects) made prior to this script, edits and submits to feat
'''


import os
import glob
import pickle
import datetime
import shutil

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

for scan in experiment['design'].keys():
    for topup in ['topUp']:#, 'noTopUp']:
        for b0 in ['noB0']:#, 'b0']:
            for HRFmodel in ['doubleGamma']:#, 'singleGamma']:
                for thr in ['2.3','2.7','3.1','3.5']:

                    nCopes = len(glob.glob(os.path.join('data/fMRI/individual/F117/210223/functional', scan, 'allRuns/topUp/noB0/doubleGamma/3.1/secondLevel.gfeat/cope*.feat')))
                    for cope in range(nCopes):
                        print(f'\n{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Scan: {scan} | topup: {topup} | b0: {b0} | HRF model: {HRFmodel} | Cope: {cope+1}/{nCopes}')

                        outDir = os.path.join('data/fMRI/group', scan, topup, b0, HRFmodel, thr, f'cope{cope+1}.gfeat')
                        os.makedirs(os.path.dirname(outDir), exist_ok=True)

                        if not os.path.isdir(outDir):

                            print('Analysis not found, analysing...')
                            designFile = 'data/fMRI/designs/thirdLevel/design.fsf'

                            with open(designFile, 'r') as file:
                                fileData = file.read()

                            # replace scan
                            fileData = fileData.replace('BFHOloc_v2', scan)

                            # replace topup and b0 types for input field and output field
                            fileData = fileData.replace('noTopUp', topup)
                            fileData = fileData.replace('noB0', b0)

                            # replace HRFmodel
                            fileData = fileData.replace('doubleGamma', HRFmodel)

                            # replace threshold
                            fileData = fileData.replace('set fmri(z_thresh) 3.1', f'set fMRI(z_thresh) {thr}')
                            fileData = fileData.replace('/3.1/', f'/{thr}/')

                            # replace cope number
                            fileData = fileData.replace('cope1', f'cope{cope+1}')

                            # check out dir in design file is correct
                            lines = fileData.splitlines()
                            for line in lines:
                                if line.startswith(f'set fmri(outputdir)'):
                                    actualOutDir = line.split()[2][1:-1]
                                    if not actualOutDir == os.path.join(os.getcwd(), outDir):
                                        print(f'targeted out dir: {os.path.join(os.getcwd(), outDir)}')
                                        print(f'actual out dir: {actualOutDir}')
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

print('Done')
