import os
import glob
import datetime
import numpy as np
import nibabel as nib
from surfer import io
from surfer import Brain

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

root = os.getcwd()
hemis = ['lh', 'rh']
fsSubjectsDir = '/mnt/HDD12TB/freesurfer/subjects'
overwriteReg = False
overwritePlot = True
ROIcopes = [8, 9, 10] # zero indexed
RGBcols = {'chartreuse': '127,255,0',
           'deepskyblue': '0,191,255',
           'orange': '255,165,0',
           'fuchsia': '255,0,255',
           'red': '255,0,0',
           'gold': '255,215,0'}

thr = 3.1
uthr = 4.5

config = {'mapType': {'zstat': {'mapPath': 'stats/zstat1',
                                'mapPathConv': 'stats/zstat1',
                                'mapPathPrep': 'stats/zstat1_thr',
                                'interp': 'trilinear'},
                      'cluster': {'mapPath': 'cluster_mask_zstat1',
                                  'mapPathConv': 'cluster_mask_zstat1',
                                  'mapPathPrep': 'cluster_mask_zstat1_bin',
                                  'interp': 'nearestneighbour'}},
            'scan': {'BFHOloc_v2': {'condNames': ['body', 'face', 'house', 'object',
                                                  'face>body', 'face>house', 'face>object', 'allConds',
                                                  'body>others', 'face>others', 'house>others', 'object>others'],
                                    'clusterCols': ['chartreuse', 'deepskyblue', 'red','gold',
                                                    'red', 'red', 'red', 'red' ,
                                                    'chartreuse', 'deepskyblue', 'orange','fuchsia']},
                     'animacyAspectRatioLoc_v2': {'condNames': ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby',
                                                                'inanimate_spiky>inanimate_stubby', 'animate>inanimate', 'spiky>stubby', 'allConds',
                                                                'animate_spiky>others', 'animate_stubby>others', 'inanimate_spiky>others', 'inanimate_stubby>others',
                                                                'animate_spiky>animate_stubby', 'NML',
                                                                'animate_stubby>inanimate_spiky','inanimate_stubby>animate_spiky','inanimate_spiky>animate_stubby', 'animate_spiky>inanimate_stubby'],
                                                  'clusterCols': ['chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                                  'red', 'red', 'red', 'red',
                                                                  'chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                                  'red','red',
                                                                  'chartreuse', 'deepskyblue', 'fuchsia', 'orange']}}}

# convert EVC mask to surface
EVCmask = '/mnt/HDD12TB/masks/Wang_2015/V1-V4.nii.gz'
for hemi in hemis:
    surfFile = f'{EVCmask[:-7]}_{hemi}.mgh'
    surfLabel = f'{EVCmask[:-7]}_{hemi}.label'
    if not os.path.isfile(surfLabel):
        freesurferCommand = f'mri_vol2surf --mov {EVCmask} --out {surfFile} --regheader fsaverage --hemi {hemi}'
        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
        freesurferCommand = f'mri_cor2label --i {surfFile} --surf fsaverage {hemi} --id 1 --l {surfLabel}'
        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
#############################################
# convert category-selective ROIs to labels #
#############################################

condNames = config['scan']['BFHOloc_v2']['condNames']
colours = config['scan']['BFHOloc_v2']['clusterCols']
runDirs = sorted(glob.glob(os.path.join('data/fMRI/group/BFHOloc_v2')))


# convert category-selective ROI activations to labels
for c in ROIcopes:

    copeDir = os.path.join('data/fMRI/group/BFHOloc_v2', topup, b0, HRFmodel, f'cope{c + 1}.gfeat/cope1.feat')
    condName = condNames[c]

    inFile = f'{copeDir}/cluster_mask_zstat1.nii.gz'
    prepFile = f'{copeDir}/cluster_mask_zstat1_highres_bin.nii.gz'

    if not os.path.exists(prepFile) or overwriteReg:
        # convert to high res
        os.system(f'fslmaths {inFile} -bin {prepFile}')  # binarize to remove different cluster labels

    for hemi in hemis:

        # convert volume to surface
        surfFile = f'{prepFile[:-7]}_{hemi}.mgh'
        if not os.path.isfile(surfFile) or overwriteReg:
            freesurferCommand = f'mri_vol2surf --mov {prepFile} --out {surfFile} --regheader {subject} --hemi {hemi}'
            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

        # convert surface to label
        surfLabel = f'{surfFile[:-4]}.label'
        if not os.path.isfile(surfLabel) or overwriteReg:
            freesurferCommand = f'mri_cor2label --i {surfFile} --surf {subject} {hemi} --id 1 --l {surfLabel}'
            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

#############################
# convert EVC mask to label #
#############################

# dilate the mask a bit for better presentation on the surface
EVCmaskDilated = 'data/fMRI/group/masks/V1-V4dilated.nii.gz'
os.makedirs(os.path.dirname(EVCmaskDilated), exist_ok=True)
os.system(f'fslmaths {EVCmask} -kernel box 3 -dilM {EVCmaskDilated}')

# convert to native highres
xmat = os.path.join(runDirs[0], 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat')
EVCmaskNative = f'data/fMRI/individual/{subject}/{session}/masks/V1-V4highres.nii.gz'
if not os.path.exists(EVCmaskNative) or overwriteReg:
    xmat = os.path.join(runDirs[0], 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat') # just use whichever xmat from the ROI transform
    os.system(f'flirt -in {EVCmask} -ref {anatomicalFile} -applyxfm -init {xmat} -out {EVCmaskNative} -interp nearestneighbour')

for hemi in hemis:

    # convert volume to surface
    surfFile = f'{EVCmaskNative[:-7]}_{hemi}.mgh'
    if not os.path.isfile(surfFile) or overwriteReg:
        freesurferCommand = f'mri_vol2surf --mov {EVCmaskNative} --out {surfFile} --regheader {subject} --hemi {hemi}'
        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

    # convert surface to label
    surfLabel = f'{EVCmaskNative[:-7]}_{hemi}.label'
    if not os.path.isfile(surfLabel) or overwriteReg:
        freesurferCommand = f'mri_cor2label --i {surfFile} --surf {subject} {hemi} --id 1 --l {surfLabel}'
        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

####################
# plot activations #
####################

for scan in experiment['design'].keys():

    condNames = config['scan'][scan]['condNames']
    colours = config['scan'][scan]['clusterCols']

    for topup in ['topUp']: #, 'noTopUp']:
        for b0 in ['noB0']: #, 'b0']:
            for HRFmodel in ['doubleGamma']: #, 'singleGamma']:

                print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Scan: {scan} | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                copeDirs = sorted(glob.glob(os.path.join('data/fMRI/group', scan, topup, b0, HRFmodel, 'cope*.gfeat')))
                nCopes = len(copeDirs)
                if nCopes != len(condNames):
                    raise Exception('Number of cope dirs does not match number of condition names')

                for mapType in config['mapType'].keys():

                    mapPath = config['mapType'][mapType]['mapPath']
                    mapPathConv = config['mapType'][mapType]['mapPathConv']
                    mapPathPrep = config['mapType'][mapType]['mapPathPrep']
                    interp = config['mapType'][mapType]['interp']

                    ### EACH CONTRAST INDIVIDUALLY
                    for c in range(nCopes):

                        copeDir = os.path.join('data/fMRI/group', scan, topup, b0, HRFmodel, f'cope{c+1}.gfeat')
                        condName = condNames[c]
                        resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, 'allSubjects', scan, condName, mapType)
                        os.makedirs(resDir, exist_ok=True)

                        inFile = f'{copeDir}/cope1.feat/{mapPath}.nii.gz'
                        prepFile = f'{copeDir}/cope1.feat/{mapPathPrep}.nii.gz'

                        if not os.path.exists(prepFile) or overwriteReg:

                            # apply further preprocessing
                            if mapType == 'zstat':
                                os.system(f'fslmaths {inFile} -thr {thr} {prepFile}')  # apply threshold to zstat images
                            if mapType == 'cluster':
                                os.system(f'fslmaths {inFile} -bin {prepFile}') # binarize cluster masks

                        for hemi in hemis:

                            # convert volume to surface
                            surfFile = f'{prepFile[:-7]}_{hemi}.mgh'
                            if not os.path.isfile(surfFile) or overwriteReg:
                                freesurferCommand = f'mri_vol2surf --mov {prepFile} --out {surfFile} --regheader fsaverage --hemi {hemi}'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                            for view in ['inferior','lateral']:
                                imageFile = os.path.join(resDir, f'{view}_{hemi}.png')
                                if not os.path.isfile(imageFile) or overwritePlot:

                                    freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/fsaverage/surf/{hemi}.inflated:curvature_method=binary'
                                    if mapType == 'zstat':
                                        freesurferCommand += f':overlay={surfFile}:overlay_threshold={thr},{uthr}'
                                    elif mapType == 'cluster':
                                        colour = RGBcols[colours[c]]
                                        freesurferCommand += f':overlay={surfFile}:overlay_custom=1,{colour}'
                                    freesurferCommand += f' -layout 1 -viewport 3d -view {view} -ss {imageFile} 1 autotrim'
                                    print(freesurferCommand)
                                    os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')



                    ### FOUR CATEGORIES OVERLAID
                    if mapType == 'cluster':
                        for contrast, copeBuffer in zip(['noContrast', 'versusOthers', 'versusOpposite'],[1,9,15]):
                            if contrast != 'versusOpposite' and scan != 'BFHOloc_v2':

                                # resolve overlap by highest z-score
                                templateData = nib.load(os.path.join(copeDirs[0], 'cope1.feat/stats/zstat1.nii.gz'))
                                volDim = list(templateData.get_fdata().shape)
                                volDim2 = volDim.copy()
                                volDim2.append(4)

                                # values to base the selection on
                                allZvals = np.zeros(volDim2)
                                for c, copeDir in enumerate(copeDirs[:4]):
                                    theseZvals = nib.load(os.path.join(copeDir, f'cope1.feat/stats/zstat1_thr.nii.gz'))
                                    allZvals[:, :, :, c] = theseZvals.get_fdata()

                                # values to be selected
                                allVals = np.zeros(volDim2)
                                for c, copeDir in enumerate(copeDirs[:4]):
                                    theseVals = nib.load(os.path.join(copeDir, f'cope1.feat/{mapPathPrep}.nii.gz'))
                                    allVals[:, :, :, c] = theseVals.get_fdata()

                                # find condition of max response for each voxel
                                maxIdx = np.argmax(allZvals, axis=3)

                                # generate new set of maps
                                newVals = np.zeros(allVals.shape)
                                for x in range(volDim[0]):
                                    for y in range(volDim[1]):
                                        for z in range(volDim[2]):
                                            thisMaxIdx = maxIdx[x, y, z]
                                            newVals[x, y, z, thisMaxIdx] = allVals[x, y, z, thisMaxIdx]

                                # save all maps in one file with value of each voxel reflecting the condition
                                allCondsFile = os.path.join(os.path.dirname(copeDir), f'customContrasts/allCondsMax_{contrast}.nii.gz')
                                os.makedirs(os.path.dirname(allCondsFile), exist_ok=True)

                                for c in range(4):

                                    # save out nifti
                                    copeDir = copeDirs[c + copeBuffer]
                                    maxMapPath = os.path.join(copeDir, f'cope1.feat/{mapPathPrep}_max.nii.gz')
                                    nib.save(nib.Nifti1Image(newVals[:, :, :, c]*(c+1), templateData.affine, templateData.header), maxMapPath)

                                    if c == 0:
                                        fslmathsCommand = f'fslmaths {maxMapPath}'
                                    else:
                                        fslmathsCommand += f' -add {maxMapPath}'

                                fslmathsCommand += f' {allCondsFile}'
                                os.system(fslmathsCommand)

                                # plot surface
                                resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel,
                                                      'allSubjects', scan, 'all', mapType, contrast)
                                os.makedirs(resDir, exist_ok=True)

                                for hemi in hemis:
                                    surfFile = f'{allCondsFile[:-7]}_{hemi}.mgh'
                                    if not os.path.isfile(surfFile) or overwritePlot:
                                        freesurferCommand = f'mri_vol2surf --mov {allCondsFile} --out {surfFile} --regheader fsaverage --hemi {hemi}'
                                        os.system(
                                            f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                    for view in ['inferior', 'lateral']:
                                        outFile = os.path.join(resDir, f'{view}_{hemi}.png')
                                        if not os.path.isfile(outFile) or overwritePlot:

                                            freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/fsaverage/surf/{hemi}.inflated' \
                                                                f':curvature_method=binary' \
                                                                f':overlay={surfFile}' \
                                                                f':overlay_custom=1,{RGBcols[colours[0]]},2,{RGBcols[colours[1]]},3,{RGBcols[colours[2]]},4,{RGBcols[colours[3]]}' \
                                                                f' -layout 1 -viewport 3d -view {view} -ss {outFile} 1 autotrim'
                                            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
