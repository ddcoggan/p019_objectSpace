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
overwriteReg = 0
overwritePlot = 0

RGBcols = {'chartreuse': '127,255,0',
           'deepskyblue': '0,191,255',
           'orange': '255,165,0',
           'fuchsia': '255,0,255',
           'red': '255,0,0',
           'gold': '255,215,0'}

thr = 3.1
uthr = 16

config = {'mapType': {'zstat': {'mapPath': 'stats/zstat1',
                                'mapPathPrep': 'stats/zstat1_thr'},
                      'cluster': {'mapPath': 'cluster_mask_zstat1',
                                  'mapPathPrep': 'cluster_mask_zstat1_bin'}},
            'scan': {'BFHOloc_v2': {'condNames': ['body', 'face', 'house', 'object',
                                                  'face>body', 'face>house', 'face>object', 'allConds',
                                                  'body>others', 'face>others', 'house>others', 'object>others'],
                                    'clusterCols': ['chartreuse', 'deepskyblue', 'red','gold',
                                                    'red', 'red', 'red', 'red' ,
                                                    'chartreuse', 'deepskyblue', 'orange','fuchsia']},
                     'animacyAspectRatioLoc_v2': {'condNames': ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby',
                                                                'inanimate_spiky>inanimate_stubby', 'animate>inanimate', 'spiky>stubby', 'allConds',
                                                                'animate_spiky>others', 'animate_stubby>others', 'inanimate_spiky>others', 'inanimate_stubby>others',
                                                                'animate_spiky>animate_stubby'],
                                                  'clusterCols': ['chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                                  'red', 'red', 'red', 'red',
                                                                  'chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                                  'red']}}}

# convert EVC mask to surface
EVCmask = '/mnt/HDD12TB/masks/Wang_2015/V1-V4.nii.gz'
EVCmasks = 'V1','V2','V3'
for EVCmask in EVCmasks:
    maskPath = f'/mnt/HDD12TB/masks/Wang_2015/{EVCmask}.nii.gz'
    for hemi in hemis:
        surfFile = f'{maskPath[:-7]}_{hemi}.mgh'
        surfLabel = f'{maskPath[:-7]}_{hemi}.label'
        if not os.path.isfile(surfLabel):
            freesurferCommand = f'mri_vol2surf --mov {maskPath} --out {surfFile} --regheader fsaverage --hemi {hemi}'
            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
            freesurferCommand = f'mri_cor2label --i {surfFile} --surf fsaverage {hemi} --id 1 --l {surfLabel}'
            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

for subject in list(experiment['scanInfo'].keys())[1:]:
    for session in experiment['scanInfo'][subject].keys():
        for scan in experiment['design'].keys():

            condNames = config['scan'][scan]['condNames']
            colours = config['scan'][scan]['clusterCols']
            runDirs = sorted(glob.glob(os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'run??')))

            for topup in ['topUp']: #, 'noTopUp']:
                for b0 in ['noB0']: #, 'b0']:
                    for HRFmodel in ['doubleGamma']: #, 'singleGamma']:

                        print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                        copeDirs = sorted(glob.glob(os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel, 'secondLevel.gfeat/cope*.feat')))
                        nCopes = len(copeDirs)

                        for mapType in config['mapType'].keys():

                            mapPath = config['mapType'][mapType]['mapPath']
                            mapPathPrep = config['mapType'][mapType]['mapPathPrep']

                            ### EACH CONTRAST INDIVIDUALLY
                            for c in range(nCopes):

                                copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel, f'secondLevel.gfeat/cope{c + 1}.feat')

                                inFile = f'{copeDir}/{mapPath}.nii.gz'
                                outFile = f'{copeDir}/{mapPathPrep}.nii.gz'
                                if not os.path.exists(outFile) or overwriteReg:
                                    if mapType == 'zstat':
                                        os.system(f'fslmaths {inFile} -thr {thr} {outFile}')  # apply threshold to zstat images
                                    if mapType == 'cluster':
                                        os.system(f'fslmaths {inFile} -bin {outFile}') # binarize cluster masks

                                condName = condNames[c]


                                resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, condName, mapType)
                                os.makedirs(resDir, exist_ok=True)
                                inFile = outFile

                                for hemi in hemis:

                                    EVClabel = f'{EVCmask[:-7]}_{hemi}.label'

                                    # convert volume to surface
                                    surfFile = f'{inFile[:-7]}_{hemi}.mgh'
                                    if not os.path.isfile(surfFile) or overwriteReg:
                                        freesurferCommand = f'mri_vol2surf --mov {inFile} --out {surfFile} --regheader fsaverage --hemi {hemi}'
                                        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                    for view in ['inferior','lateral']:
                                        outFile = os.path.join(resDir, f'{view}_{hemi}.png')
                                        if not os.path.isfile(outFile) or overwritePlot:

                                            freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/fsaverage/surf/{hemi}.inflated:curvature_method=binary'
                                            if mapType == 'zstat':
                                                freesurferCommand += f':overlay={surfFile}:overlay_threshold={thr},{uthr}'
                                            elif mapType == 'cluster':
                                                colour = RGBcols[colours[c]]
                                                freesurferCommand += f':overlay={surfFile}:overlay_custom=1,{colour}'
                                            freesurferCommand += f':label={EVClabel}:label_color=white -layout 1 -viewport 3d -view {view} -ss {outFile} 1 autotrim'
                                            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')



                            ### FOUR CATEGORIES OVERLAID
                            if mapType == 'cluster':
                                for contrast, copeBuffer in zip(['noContrast', 'versusOthers'],[1,9]):

                                    # resolve overlap by highest z-score
                                    templateData = nib.load(os.path.join(copeDirs[0], 'stats/zstat1.nii.gz'))
                                    volDim = list(templateData.get_fdata().shape)
                                    volDim2 = volDim.copy()
                                    volDim2.append(4)

                                    # values to base the selection on
                                    allZvals = np.zeros(volDim2)
                                    for c in range(4):
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + copeBuffer}.feat')
                                        theseZvals = nib.load(os.path.join(copeDir, f'stats/zstat1_thr.nii.gz'))
                                        allZvals[:, :, :, c] = theseZvals.get_fdata()

                                    # values to be selected
                                    allVals = np.zeros(volDim2)
                                    for c in range(4):
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + copeBuffer}.feat')
                                        theseVals = nib.load(os.path.join(copeDir, f'{mapPathPrep}.nii.gz'))
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
                                    allCondsFile = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/allCondsMax_{contrast}.nii.gz')

                                    for c in range(4):

                                        # save out nifti
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + copeBuffer}.feat')
                                        maxMapPath = os.path.join(copeDir, f'{mapPathPrep}_max.nii.gz')
                                        nib.save(nib.Nifti1Image(newVals[:, :, :, c]*(c+1), templateData.affine, templateData.header), maxMapPath)

                                        if c == 0:
                                            fslmathsCommand = f'fslmaths {maxMapPath}'
                                        else:
                                            fslmathsCommand += f' -add {maxMapPath}'

                                    fslmathsCommand += f' {allCondsFile}'
                                    os.system(fslmathsCommand)

                                    # plot surface
                                    resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel,
                                                          subject, scan, 'all', mapType, contrast)
                                    os.makedirs(resDir, exist_ok=True)

                                    for hemi in hemis:
                                        surfFile = f'{allCondsFile[:-7]}_{hemi}.mgh'
                                        if not os.path.isfile(surfFile) or overwritePlot:
                                            freesurferCommand = f'mri_vol2surf --mov {allCondsFile} --out {surfFile} --regheader fsaverage --hemi {hemi}'
                                            os.system(
                                                f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                        EVClabel = f'{EVCmask[:-7]}_{hemi}.label'

                                        for view in ['inferior', 'lateral']:
                                            outFile = os.path.join(resDir, f'{view}_{hemi}.png')
                                            if not os.path.isfile(outFile) or overwritePlot:

                                                freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/fsaverage/surf/{hemi}.inflated' \
                                                                    f':curvature_method=binary' \
                                                                    f':overlay={surfFile}' \
                                                                    f':overlay_custom=1,{RGBcols[colours[0]]},2,{RGBcols[colours[1]]},3,{RGBcols[colours[2]]},4,{RGBcols[colours[3]]}'
                                                freesurferCommand += f':label={EVClabel}:label_color=white -layout 1 -viewport 3d -view {view} -ss {outFile} 1 autotrim'
                                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
