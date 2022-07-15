import os
import glob
import datetime
import sys
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

# analysis 2: timecourses in V1
maskPath = '/mnt/HDD12TB/masks/Wang_2015/V1.nii.gz'
nDyn = 134
for topup in ['topUp', 'noTopUp']:
    for b0 in ['b0', 'noB0']:
        for HRFmodel in ['doubleGamma', 'singleGamma']:
            for s, scan in enumerate(experiment['design'].keys()):

                print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | V1 timecourse | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel} | Scan: {scan}')

                resDir = os.path.join('analysis/results/fMRI/dataQuality/timeCoursesV1', topup, b0, HRFmodel, scan)
                os.makedirs(resDir, exist_ok=True)
                subjects = list(experiment['scanInfo'].keys())
                data = np.empty(shape=[nDyn, len(subjects)])
                for s, subject in enumerate(experiment['scanInfo'].keys()):
                    for session in experiment['scanInfo'][subject].keys():
                        runDirs = sorted(glob.glob(os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'run??')))
                        nRuns = len(runDirs)
                        timeseries = np.empty(shape=[nDyn, nRuns])
                        for r, runDir in enumerate(runDirs):

                            # transform V1 mask
                            xmat = os.path.join(runDir, 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2example_func.mat')
                            refFile = os.path.join(runDir, 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/example_func.nii.gz')
                            outFile = os.path.join(runDir, 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/V1.nii.gz')
                            if not os.path.isfile(outFile):
                                os.system(f'flirt -in {maskPath} -ref {refFile} -applyxfm -init {xmat} -out {outFile}')

                            # get timeseries
                            timeseriesPath = os.path.join(runDir, 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/filtered_func_data.nii.gz')
                            timeseriesInfo = os.popen(f'fslmeants -i {timeseriesPath} -m {outFile}').read().split()
                            timeseries[:,r] = [float(x) for x in timeseriesInfo]

                        timeseriesMean = np.mean(timeseries, axis = 1) # average across runs
                        timeseriesMean -= np.mean(timeseriesMean) # mean center
                        data[:, s] = timeseriesMean
                df = pd.DataFrame(data, columns=subjects)
                ax = df.plot.line(title='timeseries in V1', figsize=(20,3))
                fig = ax.get_figure()
                fig.savefig(os.path.join(resDir, 'timeseries.png'))
                plt.show()
                plt.close(fig)
