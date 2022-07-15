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

# analysis 1: number of active voxels in all conditions
thr = 3.1
for topup in ['topUp', 'noTopUp']:
    for b0 in ['b0', 'noB0']:
        for HRFmodel in ['doubleGamma', 'singleGamma']:

            print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Active Voxel Count | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel}')
            resDir = os.path.join('analysis/results/fMRI/dataQuality/activeVoxels', topup, b0, HRFmodel)
            os.makedirs(resDir, exist_ok=True)
            data = []
            subjects = list(experiment['scanInfo'].keys())
            for subject in experiment['scanInfo'].keys():
                for session in experiment['scanInfo'][subject].keys():
                    voxelCounts = []
                    for s, scan in enumerate(experiment['design'].keys()):
                        zstatPath = os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'allRuns', topup, b0, HRFmodel, 'secondLevel.gfeat/cope8.feat/stats/zstat1.nii.gz')
                        nVoxels = int(os.popen(f'fslstats {zstatPath} -l {thr} -V').read().split(sep=' ')[0])
                        voxelCounts.append(nVoxels)
                    data.append(voxelCounts)
            df = pd.DataFrame(data, index=subjects, columns=['AAR', 'BFHO'])
            ax = df.plot.bar(title='number of active voxels (Z>3.1)')
            fig = ax.get_figure()
            fig.savefig(os.path.join(resDir, 'activeVoxels.png'))
            plt.show()
            plt.close(fig)

