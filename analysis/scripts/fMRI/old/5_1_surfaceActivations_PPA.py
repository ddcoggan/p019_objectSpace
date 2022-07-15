import os
import glob
import datetime
import numpy as np
import nibabel as nib
from surfer import io
from surfer import Brain

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

hemis = ['lh', 'rh']
fsSubjectsDir = '/mnt/HDD12TB/freesurfer/subjects'
overwriteReg = 0
overwritePlot = 0

RGBcols = {'orange': '255,165,0',
           'fuchsia': '255,0,255'}

colourMapMax = 50



for subject in list(experiment['scanInfo'].keys())[1:]:
    for session in experiment['scanInfo'][subject].keys():
        for topup in ['topUp']:  # , 'noTopUp']:
            for b0 in ['noB0']:  # , 'b0']:
                for HRFmodel in ['doubleGamma']:  # , 'singleGamma']:

                    resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, 'animacyAspectRatioLoc_v2', 'inanimate_spiky>inanimate_stubby')
                    os.makedirs(resDir, exist_ok=True)

                    sessDir =  f'data/fMRI/individual/{subject}/{session}'
                    houseMap = f'{sessDir}/functional/BFHOloc_v2/allRuns/{topup}/{b0}/doubleGamma/secondLevel.gfeat/cope11.feat/stats/zstat1_thr.nii.gz'
                    contrastMap = f'{sessDir}/functional/animacyAspectRatioLoc_v2/allRuns/{topup}/{b0}/doubleGamma/secondLevel.gfeat/cope5.feat/stats/cope1.nii.gz'

                    for hemi in hemis:


                        PPAest = f'{sessDir}/masks/floodfill/estimates/ppa_{hemi}.nii.gz'
                        if os.path.isfile(PPAest):

                            # create PPA mask, all voxels above thr
                            PPAmask = f'{sessDir}/masks/ppa_{hemi}_allVoxels.nii.gz'
                            if not os.path.isfile(PPAmask):
                                os.system(f'fslmaths {houseMap} -mul {PPAest} -bin {PPAmask}')

                            # convert PPA mask to surface then to label
                            PPAoverlay = f'{PPAmask[:-7]}.mgh'
                            PPAlabel = f'{PPAmask[:-7]}.label'
                            if not os.path.isfile(PPAlabel):
                                freesurferCommand = f'mri_vol2surf --mov {PPAmask} --out {PPAoverlay} --regheader fsaverage --hemi {hemi}'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                                freesurferCommand = f'mri_cor2label --i {PPAoverlay} --surf fsaverage {hemi} --id 1 --l {PPAlabel}'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                            # convert activation map to label
                            contrastOverlay= f'{contrastMap[:-7]}.mgh'
                            if not os.path.isfile(contrastOverlay):
                                freesurferCommand = f'mri_vol2surf --mov {contrastMap} --out {contrastOverlay} --regheader fsaverage --hemi {hemi}'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                            # plot
                            outFile = os.path.join(resDir, f'copeInPPA_{hemi}_inferior.png')
                            if not os.path.isfile(outFile):
                                print(f'{os.getcwd()}/{outFile}')
                                freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/fsaverage/surf/{hemi}.inflated:curvature_method=binary' \
                                                    f'overlay={contrastOverlay}:overlay_mask={PPAlabel}:overlay_threshold={-colourMapMax},{colourMapMax}' \
                                                    f':overlay_custom={-colourMapMax},{RGBcols["fuchsia"]},{colourMapMax},{RGBcols["orange"]} ' \
                                                    f'-layout 1 -viewport 3d -view inferior -ss {outFile} 1 autotrim --noquit'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')