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
thr = 3.1
uthr = 16
thrBin = 0.5
uthrBin = 2
hemis = ['lh', 'rh']
distance = 800
fsSubjectsDir = '/mnt/HDD12TB/freesurfer/subjects'
overwriteReg = 0
overwritePlot = 0

RGBcols = {'chartreuse': '127,255,0',
           'deepskyblue': '0,191,255',
           'orange': '255,165,0',
           'fuchsia': '255,0,255',
           'red': '255,0,0',
           'gold': '255,215,0'}

config = {'BFHOloc_v2': {'condNames': ['body', 'face', 'house', 'object',
                                       'face>body', 'face>house', 'face>object', 'allConds',
                                       'body>others', 'face>others', 'house>others', 'object>others'],
                         'clusterCols': ['chartreuse', 'deepskyblue', 'red','gold',
                                         'red', 'red', 'red', 'red' ,
                                         'chartreuse', 'deepskyblue', 'orange','fuchsia']},
          'animacyAspectRatioLoc_v2': {'condNames': ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby',
                                                     'inanimate_spiky>inanimate_stubby', 'animate>inanimate', 'spiky>stubby', 'allConds',
                                                     'animate_spiky>others', 'animate_stubby>others', 'inanimate_spiky>others', 'inanimate_stubby>others',
                                                     'animate_spiky>animate_stubby','NML','animate_spiky>inanimate_stubby',
                                                     'animate_stubby>inanimate_spiky','inanimate_stubby>animate_spiky','inanimate_spiky>animate_stubby'],
                                       'clusterCols': ['chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                       'red', 'red', 'red', 'red',
                                                       'chartreuse', 'deepskyblue', 'orange', 'fuchsia',
                                                       'red', 'red',
                                                       'chartreuse', 'deepskyblue', 'fuchsia', 'orange']}}

for subject in list(experiment['scanInfo'].keys())[1:]:
    anatomicalFile = f'/mnt/HDD12TB/freesurfer/subjects/{subject}/mri/anatomical.nii'
    for session in experiment['scanInfo'][subject].keys():
        for scan in experiment['design'].keys():
            runDirs = sorted(glob.glob(os.path.join('data/fMRI/individual', subject, session, 'functional', scan, 'run??')))
            condNames = config[scan]['condNames']
            colours = config[scan]['clusterCols']
            for topup in ['topUp']:
                for b0 in ['noB0']:
                    for HRFmodel in ['doubleGamma']:

                        print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | topUp: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                        copeDirs = sorted(glob.glob(os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,'secondLevel.gfeat/cope*.feat')))
                        nCopes = len(copeDirs)
                        for c in range(nCopes):

                            copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel, f'secondLevel.gfeat/cope{c+1}.feat')

                            # convert zstat to high res
                            inFile = os.path.join(copeDir, 'stats/zstat1.nii.gz')
                            xmat = os.path.join(runDirs[0], 'outputs', topup, b0, HRFmodel, 'firstLevel.feat/reg/standard2highres.mat')
                            outFile = os.path.join(copeDir, 'stats/zstat1_highres.nii.gz')

                            if not os.path.exists(outFile) or overwritePlot:
                                os.system(f'flirt -in {inFile} -ref {anatomicalFile} -applyxfm -init {xmat} -out {outFile}')
                                os.system(f'fslmaths {outFile} -thr {thr} {outFile}')

                            # convert high res volume to surface
                            inFile = outFile
                            condName = condNames[c]
                            resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, condName)
                            if not os.path.isdir(resDir):
                                os.makedirs(resDir, exist_ok=True)

                            for hemi in hemis:
                                surfFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}.mgh')
                                if not os.path.isfile(outFile) or overwriteReg:
                                    freesurferCommand = f'mri_vol2surf --mov {inFile} --out {surfFile} --regheader {subject} --hemi {hemi}'
                                    os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                for view in ['inferior', 'lateral']:

                                    imageFile = os.path.join(resDir, f'{view}_{hemi}.png')

                                    if not os.path.isfile(imageFile) or overwritePlot:

                                        # plot activations on surface (each contrast separately)
                                        freesurferCommand = f'freeview -f /mnt/HDD12TB/freesurfer/subjects/{subject}/surf/{hemi}.inflated:curvature_method=binary'
                                        colour = RGBcols[colours[c]]
                                        freesurferCommand += f':overlay={surfFile}:overlay_custom=1,255,0,0,16,255,255,0: -layout 1 -viewport 3d -view {view} 1 autotrim'
                                        os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
'''
                        ### EACH CATEGORY RAW PARAMETER ESTIMATES
                        # plot activations on surface (each category raw beta values)
                        resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'PEs', 'overlapping')
                        for hemi in hemis:
                            finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                            if not os.path.exists(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 1}.feat')
                                    condName = condNames[c]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}.mgh'))
                                    brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c]
                                os.makedirs(resDir, exist_ok=True)
                                brain.show_view(distance=distance, view='ventral')
                                brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                brain.show_view(distance=distance, view='lateral')
                                brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                        # plot which contrasts have max value
                        copeImg = nib.load(os.path.join(copeDirs[0], 'stats/zstat1_highres.nii.gz'))
                        volDim = list(copeImg.get_fdata().shape)
                        volDimCope = volDim.copy()
                        volDimCope.append(4)
                        allCopes = np.zeros(volDimCope)
                        for c in range(4):
                            copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                   f'secondLevel.gfeat/cope{c + 1}.feat')
                            copeImg = nib.load(os.path.join(copeDir, 'stats/zstat1_highres.nii.gz'))
                            allCopes[:,:,:,c] = copeImg.get_fdata()
                        maxIdx = np.argmax(allCopes, axis = 3)
                        maxCopes = np.zeros(allCopes.shape)
                        for x in range(volDim[0]):
                            for y in range(volDim[1]):
                                for z in range(volDim[2]):
                                    maxCope = maxIdx[x,y,z]
                                    maxCopes[x,y,z,maxCope] = allCopes[x,y,z,maxCope]
                        minCopeVals = np.amin(allCopes, axis = 3)
                        for c in range(4):
                            copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                   f'secondLevel.gfeat/cope{c + 1}.feat')
                            maxCopeFile = os.path.join(copeDir, 'stats/zstat1_highres_max.nii.gz')
                            nib.save(nib.Nifti1Image(maxCopes[:,:,:,c], copeImg.affine, copeImg.header), maxCopeFile)
                            maxCopeFileBin = os.path.join(copeDir, 'stats/zstat1_highres_max_bin.nii.gz')
                            os.system(f'fslmaths {maxCopeFile} -bin {maxCopeFileBin}')
                            condName = condNames[c]
                            for hemi in hemis:
                                # retaining PEs
                                outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh')
                                if not os.path.isfile(outFile) or overwritePlot:
                                    freesurferCommand = f'mri_vol2surf --mov {maxCopeFile} --out {outFile} --regheader {subject} --hemi {hemi}'
                                    os.system(
                                        f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                                # PE map binarized
                                outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh')
                                if not os.path.isfile(outFile) or overwritePlot:
                                    freesurferCommand = f'mri_vol2surf --mov {maxCopeFileBin} --out {outFile} --regheader {subject} --hemi {hemi}'
                                    os.system(
                                        f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                        resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'PEs', 'nonoverlapping')
                        for hemi in hemis:
                            # retaining PEs
                            finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                            if not os.path.isfile(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",
                                              subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 1}.feat')
                                    condName = condNames[c]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh'))
                                    
                                    brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c]
                                os.makedirs(resDir, exist_ok=True)
                                brain.show_view(distance=distance, view='ventral')
                                brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                brain.show_view(distance=distance, view='lateral')
                                brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                            # binarize PEs
                            finalFile = os.path.join(resDir, f'lateral_{hemi}_bin.png')
                            if not os.path.isfile(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",
                                              subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 1}.feat')
                                    condName = condNames[c]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh'))
                                    brain.add_overlay(zstatSurf, thrBin, uthrBin, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c]
                                os.makedirs(resDir, exist_ok=True)
                                brain.show_view(distance=distance, view='ventral')
                                brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}_bin.png'))
                                brain.show_view(distance=distance, view='lateral')
                                brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}_bin.png'))

                        ### EACH CATEGORY VERSUS OTHERS
                        # plot activations on surface
                        resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'versus_others', 'overlapping')
                        for hemi in hemis:
                            finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                            if not os.path.isfile(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 9}.feat')
                                    condName = condNames[c+8]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}.mgh'))
                                    
                                    brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c+8]
                                    os.makedirs(resDir, exist_ok=True)
                                    brain.show_view(distance=distance, view='ventral')
                                    brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                    brain.show_view(distance=distance, view='lateral')
                                    brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                        # plot which contrasts have max value
                        copeImg = nib.load(os.path.join(copeDirs[0], 'stats/zstat1_highres.nii.gz'))
                        volDim = list(copeImg.get_fdata().shape)
                        volDimCope = volDim.copy()
                        volDimCope.append(4)
                        allCopes = np.zeros(volDimCope)
                        for c in range(4):
                            copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                   f'secondLevel.gfeat/cope{c + 9}.feat')
                            copeImg = nib.load(os.path.join(copeDir, 'stats/zstat1_highres.nii.gz'))
                            allCopes[:,:,:,c] = copeImg.get_fdata()
                        maxIdx = np.argmax(allCopes, axis = 3)
                        maxCopes = np.zeros(allCopes.shape)
                        for x in range(volDim[0]):
                            for y in range(volDim[1]):
                                for z in range(volDim[2]):
                                    maxCope = maxIdx[x,y,z]
                                    maxCopes[x,y,z,maxCope] = allCopes[x,y,z,maxCope]
                        minCopeVals = np.amin(allCopes, axis = 3)
                        for c in range(4):
                            copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                   f'secondLevel.gfeat/cope{c + 9}.feat')
                            maxCopeFile = os.path.join(copeDir, 'stats/zstat1_highres_max.nii.gz')
                            nib.save(nib.Nifti1Image(maxCopes[:,:,:,c], copeImg.affine, copeImg.header), maxCopeFile)
                            maxCopeFileBin = os.path.join(copeDir, 'stats/zstat1_highres_max_bin.nii.gz')
                            os.system(f'fslmaths {maxCopeFile} -bin {maxCopeFileBin}')
                            condName = condNames[c+8]
                            for hemi in hemis:
                                # retaining PEs
                                outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh')
                                if not os.path.isfile(outFile) or overwritePlot:
                                    freesurferCommand = f'mri_vol2surf --mov {maxCopeFile} --out {outFile} --regheader {subject} --hemi {hemi}'
                                    os.system(
                                        f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                                # PE map binarized
                                outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh')
                                if not os.path.isfile(outFile) or overwritePlot:
                                    freesurferCommand = f'mri_vol2surf --mov {maxCopeFileBin} --out {outFile} --regheader {subject} --hemi {hemi}'
                                    os.system(
                                        f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                        resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'versus_others', 'nonoverlapping')
                        for hemi in hemis:

                            # retaining PEs
                            finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                            if not os.path.isfile(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",
                                              subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 9}.feat')
                                    condName = condNames[c+8]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh'))
                                    
                                    brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c+8]
                                os.makedirs(resDir, exist_ok=True)
                                brain.show_view(distance=distance, view='ventral')
                                brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                brain.show_view(distance=distance, view='lateral')
                                brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                            # binarised PEs
                            finalFile = os.path.join(resDir, f'lateral_{hemi}_bin.png')
                            if not os.path.isfile(finalFile) or overwritePlot:
                                brain = Brain(subject, hemi, "inflated", background="white",
                                              subjects_dir=fsSubjectsDir)
                                for c in range(4):
                                    copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                           f'secondLevel.gfeat/cope{c + 9}.feat')
                                    condName = condNames[c+8]
                                    zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh'))
                                    brain.add_overlay(zstatSurf, thrBin, uthrBin, name=condName, sign='pos')
                                    brain.overlays[condName].pos_bar.lut_mode = colours[c+8]
                                os.makedirs(resDir, exist_ok=True)
                                brain.show_view(distance=distance, view='ventral')
                                brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}_bin.png'))
                                brain.show_view(distance=distance, view='lateral')
                                brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}_bin.png'))

                        ### EACH CATEGORY VERSUS OPPOSITE
                        if scan == 'animacyAspectRatioLoc_v2':
                            # plot activations on surface
                            resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'versus_opposite', 'overlapping')
                            for hemi in hemis:
                                finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                                if not os.path.isfile(finalFile) or overwritePlot:
                                    brain = Brain(subject, hemi, "inflated", background="white", subjects_dir=fsSubjectsDir)
                                    for c in range(4):
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + 15}.feat')
                                        condName = condNames[c + 14]
                                        zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}.mgh'))

                                        brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                        brain.overlays[condName].pos_bar.lut_mode = colours[c + 14]
                                        os.makedirs(resDir, exist_ok=True)
                                        brain.show_view(distance=distance, view='ventral')
                                        brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                        brain.show_view(distance=distance, view='lateral')
                                        brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                            # plot which contrasts have max value
                            copeImg = nib.load(os.path.join(copeDirs[0], 'stats/zstat1_highres.nii.gz'))
                            volDim = list(copeImg.get_fdata().shape)
                            volDimCope = volDim.copy()
                            volDimCope.append(4)
                            allCopes = np.zeros(volDimCope)
                            for c in range(4):
                                copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                       f'secondLevel.gfeat/cope{c + 15}.feat')
                                copeImg = nib.load(os.path.join(copeDir, 'stats/zstat1_highres.nii.gz'))
                                allCopes[:, :, :, c] = copeImg.get_fdata()
                            maxIdx = np.argmax(allCopes, axis=3)
                            maxCopes = np.zeros(allCopes.shape)
                            for x in range(volDim[0]):
                                for y in range(volDim[1]):
                                    for z in range(volDim[2]):
                                        maxCope = maxIdx[x, y, z]
                                        maxCopes[x, y, z, maxCope] = allCopes[x, y, z, maxCope]
                            minCopeVals = np.amin(allCopes, axis=3)
                            for c in range(4):
                                copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                       f'secondLevel.gfeat/cope{c + 15}.feat')
                                maxCopeFile = os.path.join(copeDir, 'stats/zstat1_highres_max.nii.gz')
                                nib.save(nib.Nifti1Image(maxCopes[:, :, :, c], copeImg.affine, copeImg.header), maxCopeFile)
                                maxCopeFileBin = os.path.join(copeDir, 'stats/zstat1_highres_max_bin.nii.gz')
                                os.system(f'fslmaths {maxCopeFile} -bin {maxCopeFileBin}')
                                condName = condNames[c + 14]
                                for hemi in hemis:
                                    # retaining PEs
                                    outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh')
                                    if not os.path.isfile(outFile) or overwritePlot:
                                        freesurferCommand = f'mri_vol2surf --mov {maxCopeFile} --out {outFile} --regheader {subject} --hemi {hemi}'
                                        os.system(
                                            f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                                    # PE map binarized
                                    outFile = os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh')
                                    if not os.path.isfile(outFile) or overwritePlot:
                                        freesurferCommand = f'mri_vol2surf --mov {maxCopeFileBin} --out {outFile} --regheader {subject} --hemi {hemi}'
                                        os.system(
                                            f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                            resDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', 'versus_opposite', 'nonoverlapping')
                            for hemi in hemis:

                                # retaining PEs
                                finalFile = os.path.join(resDir, f'lateral_{hemi}.png')
                                if not os.path.isfile(finalFile) or overwritePlot:
                                    brain = Brain(subject, hemi, "inflated", background="white",
                                                  subjects_dir=fsSubjectsDir)
                                    for c in range(4):
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + 15}.feat')
                                        condName = condNames[c + 14]
                                        zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max.mgh'))

                                        brain.add_overlay(zstatSurf, thr, uthr, name=condName, sign='pos')
                                        brain.overlays[condName].pos_bar.lut_mode = colours[c + 14]
                                    os.makedirs(resDir, exist_ok=True)
                                    brain.show_view(distance=distance, view='ventral')
                                    brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}.png'))
                                    brain.show_view(distance=distance, view='lateral')
                                    brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}.png'))

                                # binarised PEs
                                finalFile = os.path.join(resDir, f'lateral_{hemi}_bin.png')
                                if not os.path.isfile(finalFile) or overwritePlot:
                                    brain = Brain(subject, hemi, "inflated", background="white",
                                                  subjects_dir=fsSubjectsDir)
                                    for c in range(4):
                                        copeDir = os.path.join(os.path.dirname(runDirs[0]), 'allRuns', topup, b0, HRFmodel,
                                                               f'secondLevel.gfeat/cope{c + 15}.feat')
                                        condName = condNames[c + 14]
                                        zstatSurf = io.read_scalar_data(os.path.join(copeDir, f'stats/zstat1_highres_{hemi}_max_bin.mgh'))
                                        brain.add_overlay(zstatSurf, thrBin, uthrBin, name=condName, sign='pos')
                                        brain.overlays[condName].pos_bar.lut_mode = colours[c + 14]
                                    os.makedirs(resDir, exist_ok=True)
                                    brain.show_view(distance=distance, view='ventral')
                                    brain.save_image(filename=os.path.join(resDir, f'ventral_{hemi}_bin.png'))
                                    brain.show_view(distance=distance, view='lateral')
                                    brain.save_image(filename=os.path.join(resDir, f'lateral_{hemi}_bin.png'))

'''
