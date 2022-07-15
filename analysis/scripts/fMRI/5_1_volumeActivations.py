import os
import glob
import datetime
import numpy as np
import nibabel as nib
from surfer import io
from surfer import Brain
import time
#time.sleep(7200)
# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

root = os.getcwd()
hemis = ['lh', 'rh']
fsSubjectsDir = '/mnt/HDD12TB/freesurfer/subjects'
overwriteReg = True
overwritePlot = True
ROIcopes = [8, 9, 10]
ROIcols = ['yellow', 'blue', 'red']
RGBcols = {'chartreuse': '127,255,0',
           'deepskyblue': '0,191,255',
           'orange': '255,165,0',
           'fuchsia': '255,0,255',
           'red': '255,0,0',
           'blue': '0,0,255',
           'gold': '255,215,0'}
EVCmask = '/mnt/HDD12TB/masks/Wang_2015/V1-V4.nii.gz'
EVCcol = 'black'
thrs = [2.3,2.7,3.1,3.5]
uthr = 10

config = {'BFHOloc_v2': {'condNames': ['body', 'face', 'house', 'object',
                                     'face>body', 'face>house', 'face>object', 'allConds',
                                     'body>others', 'face>others', 'house>others', 'object>others'],
                       'clusterCols': ['gold', 'blue', 'red', 'chartreuse',
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
                                                     'chartreuse', 'deepskyblue', 'fuchsia', 'orange']}} # last two conditions are reversed, see fsl design


def getCopeDir(scanDir, subject, topup, b0, HRFmodel, thr, c):
    if subject != 'fsaverage':
        if c == '*':
            return glob.glob(f'{scanDir}/allRuns/{topup}/{b0}/{HRFmodel}/{str(thr)}/secondLevel.gfeat/cope*.feat')
        else:
            return f'{scanDir}/allRuns/{topup}/{b0}/{HRFmodel}/{str(thr)}/secondLevel.gfeat/cope{c + 1}.feat'
    else:
        if c == '*':
            return glob.glob(f'{scanDir}/{topup}/{b0}/{HRFmodel}/{str(thr)}/cope*.gfeat/cope1.feat')
        else:
            return f'{scanDir}/{topup}/{b0}/{HRFmodel}/{str(thr)}/cope{c + 1}.gfeat/cope1.feat'


for thr in thrs:
    for topup in ['topUp']:  # , 'noTopUp']:
        for b0 in ['noB0']:  # , 'b0']:
            for HRFmodel in ['doubleGamma']:  # , 'singleGamma']:
                subjects = list(experiment['scanInfo'].keys())[1:] + ['fsaverage']
                for subject in subjects:

                    # individual and group mean require different processing
                    if subject != 'fsaverage':
                        subjectString = subject
                        session = list(experiment['scanInfo'][subject].keys())[0]
                        preScanDir = f'data/fMRI/individual/{subject}/{session}/functional'
                        fileConfig = {'zstat': {'mapPath': 'stats/zstat1',
                                                'mapPathMid': 'stats/zstat1_highres',
                                                'mapPathFinal': 'stats/zstat1_highres_thr'},
                                      'cluster': {'mapPath': 'cluster_mask_zstat1',
                                                  'mapPathMid': 'cluster_mask_zstat1_bin',
                                                  'mapPathFinal': 'cluster_mask_zstat1_bin_highres'}}
                        anatomicalFile = f'data/fMRI/individual/{subject}/{session}/anatomical/anatomical.nii'
                    else:
                        subjectString = 'allSubjects'
                        preScanDir = 'data/fMRI/group/'
                        fileConfig = {'zstat': {'mapPath': 'stats/zstat1',
                                                'mapPathFinal': 'stats/zstat1_thr'},
                                      'cluster': {'mapPath': 'cluster_mask_zstat1',
                                                  'mapPathMid': 'cluster_mask_zstat1_bin',
                                                  'mapPathFinal': 'cluster_mask_zstat1_bin'}}




                    #############################################
                    # convert category-selective ROIs to labels #
                    #############################################

                    condNames = config['BFHOloc_v2']['condNames']
                    scanDir = preScanDir + '/BFHOloc_v2'

                    # convert category-selective ROI activations to labels
                    for c in ROIcopes:

                        condName = condNames[c]
                        copeDir = getCopeDir(scanDir, subject, topup, b0, HRFmodel, thr, c)
                        inFile = f'{copeDir}/{fileConfig["cluster"]["mapPath"]}.nii.gz'

                        # binarise to  remove different cluster labels
                        outFile = f'{copeDir}/{fileConfig["cluster"]["mapPathMid"]}.nii.gz'
                        os.system(f'fslmaths {inFile} -bin {outFile}')
                        inFile = outFile

                        # individual data converted to high res, group mean left in MNI space
                        if subject != 'fsaverage':

                            outFile = f'{copeDir}/{fileConfig["cluster"]["mapPathFinal"]}.nii.gz'
                            if not os.path.exists(outFile) or overwriteReg:

                                # convert to high res
                                xmat = os.path.join(scanDir, 'run01/outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat')
                                os.system(f'flirt -in {inFile} -ref {anatomicalFile} -applyxfm -init {xmat} -out {outFile} -interp trilinear')
                                # threshold and binarise (trilinear then thr then bin is smoother than nearest neighbour)
                                os.system(f'fslmaths {outFile} -thr 0.5 -bin {outFile}')

                                inFile = outFile

                        for hemi in hemis:

                            # convert volume to surface
                            surfFile = f'{inFile[:-7]}_{hemi}.mgh'
                            if not os.path.isfile(surfFile) or overwriteReg:
                                freesurferCommand = f'mri_vol2surf --mov {inFile} --out {surfFile} --regheader {subject} --hemi {hemi}'
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
                    os.system(f'fslmaths {EVCmask} -kernel box 4 -dilM {EVCmaskDilated}')

                    # convert to native highres unless group data
                    if subject != 'fsaverage':
                        EVCmaskNative = f'data/fMRI/individual/{subject}/{session}/masks/V1-V4highres.nii.gz'
                        if not os.path.exists(EVCmaskNative) or overwriteReg:
                            xmat = os.path.join(scanDir, 'run01/outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat') # just use whichever xmat from the ROI transform
                            os.system(f'flirt -in {EVCmaskDilated} -ref {anatomicalFile} -applyxfm -init {xmat} -out {EVCmaskNative} -interp nearestneighbour')
                    else:
                        EVCmaskNative = EVCmaskDilated

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

                        print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Plotting Activations | Subject: {subject} | Session: {session} | Scan: {scan} | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel} | thr: {thr}')

                        condNames = config[scan]['condNames']
                        colours = config[scan]['clusterCols']

                        scanDir = os.path.join(preScanDir, scan)
                        nCopes = len(condNames)

                        for mapType in fileConfig:

                            ### EACH CONTRAST INDIVIDUALLY
                            for c in range(nCopes):

                                copeDir = getCopeDir(scanDir, subject, topup, b0, HRFmodel, str(thr), c)
                                condName = condNames[c]


                                finalFile = f'{copeDir}/{fileConfig[mapType]["mapPathFinal"]}.nii.gz'

                                if not os.path.exists(finalFile) or overwriteReg:


                                    # zstat processing
                                    if mapType == 'zstat':
                                        inFile = f'{copeDir}/{fileConfig[mapType]["mapPath"]}.nii.gz'

                                        # convert to high res (ind maps only)
                                        if subject != 'fsaverage':
                                            outFile = f'{copeDir}/{fileConfig[mapType]["mapPathMid"]}.nii.gz'
                                            xmat = os.path.join(scanDir, 'run01/outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat')
                                            os.system(f'flirt -in {inFile} -ref {anatomicalFile} -applyxfm -init {xmat} -out {outFile} -interp trilinear')
                                            inFile = outFile

                                        # threshold
                                        outFile = f'{copeDir}/{fileConfig[mapType]["mapPathFinal"]}.nii.gz'
                                        os.system(f'fslmaths {inFile} -thr {thr} {outFile}')  # apply threshold to zstat images

                                    elif mapType == 'cluster':
                                        inFile = f'{copeDir}/{fileConfig[mapType]["mapPath"]}.nii.gz'
                                        outFile = f'{copeDir}/{fileConfig[mapType]["mapPathMid"]}.nii.gz'

                                        # binarise cluster masks
                                        os.system(f'fslmaths {inFile} -bin {outFile}')  # this removes the cluster labels
                                        inFile = outFile

                                        # convert to high res (ind maps only)
                                        outFile = f'{copeDir}/{fileConfig[mapType]["mapPathFinal"]}.nii.gz'
                                        if subject != 'fsaverage':
                                            xmat = os.path.join(scanDir, 'run01/outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat')
                                            os.system(f'flirt -in {inFile} -ref {anatomicalFile} -applyxfm -init {xmat} -out {outFile} -interp trilinear')
                                            os.system(f'fslmaths {outFile} -thr 0.5 -bin {outFile}')  # binarize cluster masks (trilinear then thr then bin is smoother than nearest neighbour)

                                # plot each hemisphere
                                for hemi in hemis:

                                    # get ROI labels
                                    ROIlabelPaths = []

                                    mapPathFinalROI = fileConfig['cluster']['mapPathFinal']
                                    for x in ROIcopes:
                                        scanDirROI = os.path.join(preScanDir, 'BFHOloc_v2')
                                        copeDirROI = getCopeDir(scanDirROI, subject, topup, b0, HRFmodel, str(thr), x)
                                        ROIlabelPaths.append(f'{copeDirROI}/{mapPathFinalROI}_{hemi}.label')

                                    EVClabelPath = f'{EVCmaskNative[:-7]}_{hemi}.label'

                                    # convert volume to surface
                                    surfFile = f'{inFile[:-7]}_{hemi}.mgh'
                                    if not os.path.isfile(surfFile) or overwriteReg:

                                        freesurferCommand = f'mri_vol2surf --mov {inFile} --out {surfFile} --regheader {subject} --hemi {hemi}'
                                        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                    for view in ['inferior','lateral']:
                                        for ROIsOnOff,ROIfileString in zip([False,True],['withoutROIs','withROIs']):

                                            resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subjectString, scan, condName, mapType, ROIfileString, str(thr))
                                            os.makedirs(resDir, exist_ok=True)
                                            imageFile = os.path.join(resDir, f'{view}_{hemi}.png')

                                            if not os.path.isfile(imageFile) or overwritePlot:

                                                freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/{subject}/surf/{hemi}.inflated:curvature_method=binary'

                                                if mapType == 'zstat':
                                                    freesurferCommand += f':overlay={surfFile}:overlay_custom={thr},0,255,0,{uthr},255,255,255'
                                                elif mapType == 'cluster':
                                                    colour = RGBcols[colours[c]]
                                                    freesurferCommand += f':overlay={surfFile}:overlay_custom=1,{colour}'

                                                if ROIsOnOff == True:
                                                    for ROIlabelPath, rcol in zip(ROIlabelPaths, ROIcols):
                                                        if os.path.isfile(ROIlabelPath):
                                                            freesurferCommand += f':label={ROIlabelPath}:label_color={rcol}:label_outline=yes'
                                                    freesurferCommand += f':label={EVClabelPath}:label_color={EVCcol}:label_outline=yes'

                                                freesurferCommand += f' -layout 1 -viewport 3d -view {view} -ss {imageFile} 1 autotrim'

                                                if condName in ['NML','spiky>stubby']:
                                                    os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')


                            ### FOUR CATEGORIES OVERLAID
                            if mapType == 'cluster':
                                for contrast, copeBuffer in zip(['noContrast', 'versusOthers', 'versusOpposite'],[0,8,14]):
                                    if contrast != 'versusOpposite' or scan != 'BFHOloc_v2':

                                        # resolve overlap by highest z-score
                                        if subject != 'fsaverage':
                                            templateData = nib.load(os.path.join(copeDir, 'stats/zstat1_highres.nii.gz'))
                                        else:
                                            templateData = nib.load('/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz')

                                        volDim = list(templateData.get_fdata().shape)
                                        volDim2 = volDim.copy()
                                        volDim2.append(4)

                                        # z values to base the selection on
                                        allZvals = np.zeros(volDim2)
                                        for c in range(4):
                                            copeDir = getCopeDir(scanDir, subject, topup, b0, HRFmodel, str(thr), c + copeBuffer)
                                            theseZvals = nib.load(os.path.join(copeDir, f'{fileConfig["zstat"]["mapPathFinal"]}.nii.gz'))
                                            allZvals[:, :, :, c] = theseZvals.get_fdata()

                                        # cluster values to be selected
                                        allVals = np.zeros(volDim2)
                                        for c in range(4):
                                            copeDir = getCopeDir(scanDir, subject, topup, b0, HRFmodel, str(thr), c + copeBuffer)
                                            theseVals = nib.load(os.path.join(copeDir, f'{fileConfig["cluster"]["mapPathFinal"]}.nii.gz'))
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
                                        allCondsFile = os.path.join(os.path.dirname(copeDir), f'allCondsMax_{contrast}.nii.gz')

                                        for c in range(4):

                                            # save out nifti
                                            copeDir = getCopeDir(scanDir, subject, topup, b0, HRFmodel, str(thr), c + copeBuffer)
                                            maxMapPath = f'{copeDir}/{fileConfig[mapType]["mapPathFinal"]}_max.nii.gz'
                                            nib.save(nib.Nifti1Image(newVals[:, :, :, c]*(c+1), templateData.affine, templateData.header), maxMapPath)

                                            if c == 0:
                                                fslmathsCommand = f'fslmaths {maxMapPath}'
                                            else:
                                                fslmathsCommand += f' -add {maxMapPath}'

                                        fslmathsCommand += f' {allCondsFile}'
                                        os.system(fslmathsCommand)

                                        # plot surface
                                        for hemi in hemis:

                                            # get ROI labels
                                            ROIlabelPaths = []

                                            mapPathPrepROI = fileConfig['cluster']['mapPathFinal']
                                            for x in ROIcopes:
                                                scanDirROI = os.path.join(preScanDir, 'BFHOloc_v2')
                                                copeDirROI = getCopeDir(scanDirROI, subject, topup, b0, HRFmodel, str(thr), x)
                                                ROIlabelPaths.append(f'{copeDirROI}/{mapPathPrepROI}_{hemi}.label')

                                            EVClabelPath = f'{EVCmaskNative[:-7]}_{hemi}.label'

                                            surfFile = f'{allCondsFile[:-7]}_{hemi}.mgh'
                                            if not os.path.isfile(surfFile) or overwritePlot:
                                                freesurferCommand = f'mri_vol2surf --mov {allCondsFile} --out {surfFile} --regheader {subject} --hemi {hemi}'
                                                os.system(
                                                    f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                            for view in ['inferior', 'lateral']:
                                                for ROIsOnOff, ROIfileString in zip([False, True], ['withoutROIs', 'withROIs']):
                                                    resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel,
                                                                          subjectString, scan, 'all', mapType, contrast, ROIfileString, str(thr))
                                                    os.makedirs(resDir, exist_ok=True)
                                                    imageFile = os.path.join(resDir, f'{view}_{hemi}.png')
                                                    if not os.path.isfile(imageFile) or overwritePlot:

                                                        freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/{subject}/surf/{hemi}.inflated' \
                                                                            f':curvature_method=binary' \
                                                                            f':overlay={surfFile}' \
                                                                            f':overlay_custom=1,{RGBcols[colours[copeBuffer]]},2,{RGBcols[colours[copeBuffer+1]]},3,{RGBcols[colours[copeBuffer+2]]},4,{RGBcols[colours[copeBuffer+3]]}'

                                                        if ROIsOnOff == True:
                                                            for ROIlabelPath, rcol in zip(ROIlabelPaths, ROIcols):
                                                                if os.path.isfile(ROIlabelPath):
                                                                    freesurferCommand += f':label={ROIlabelPath}:label_color={rcol}:label_outline=yes'
                                                            freesurferCommand += f':label={EVClabelPath}:label_color={EVCcol}:label_outline=yes'

                                                        freesurferCommand += f' -layout 1 -viewport 3d -view {view} -ss {imageFile} 1 autotrim'

                                                        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
