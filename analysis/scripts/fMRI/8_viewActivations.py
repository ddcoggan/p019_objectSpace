'''
loads up fsleyes and places t maps on surface so ROIs can be drawn
'''
import os

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment
subjects = list(experiment['scanInfo'].keys())
subject = subjects[2]
subject = 'F111'
print(subject)
session = 0
thr = 3.1
uthr = 6
sess = list(experiment['scanInfo'][subject].keys())[session]
sessDir = f'data/fMRI/individual/{subject}/{sess}'
maskDir = f'{sessDir}/masks/floodfill/estimates'
os.makedirs(maskDir, exist_ok=True)

# for animacy/aspect ratio
expDir = f'{sessDir}/functional/animacyAspectRatioLoc_v2/allRuns/topUp/noB0/doubleGamma/3.1/secondLevel.gfeat'
locDir = f'{sessDir}/functional/BFHOloc_v2/allRuns/topUp/noB0/doubleGamma/3.1/secondLevel.gfeat'
exampleFunc = f'{locDir}/cope1.feat/example_func.nii.gz'
evcPath = '/mnt/HDD12TB/masks/Wang_2015/V1-V4.nii.gz'
nmlPath = f'{expDir}/cope14.feat/stats/tstat1.nii.gz'
bodyPath = f'{locDir}/cope9.feat/stats/tstat1.nii.gz'
facePath = f'{locDir}/cope10.feat/stats/tstat1.nii.gz'
housePath = f'{locDir}/cope11.feat/stats/tstat1.nii.gz'
os.system(f'fsleyes  {exampleFunc} '
          f'{nmlPath} -dr {thr} {uthr} -cm green '
          f'{evcPath} -dr 0 1 -cm hot -a 30 '
          f'{bodyPath} -dr {thr} {thr} -cm yellow -a 30 '
          f'{facePath} -dr {thr} {thr} -cm blue -a 30 '
          f'{housePath} -dr {thr} {thr} -cm red -a 30')