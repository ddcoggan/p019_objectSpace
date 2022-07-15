import numpy as np
import nibabel as nib
from scipy.stats import sem
import os
import glob
import matplotlib.pyplot as plt
import datetime
import shutil
import pandas as pd
import pickle as pkl

overwrite = 0
MRthresh = 50000
regScan = 'BFHOloc_v2' # directory from which to obtain registration files
regRun = 'run05'

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

# store collated data in a dictionary of lists
maxSessions = 1
subjects = list(experiment['scanInfo'].keys())[1:]
nSubjects = len(subjects)
dqDir = 'analysis/results/fMRI/dataQuality'
os.makedirs(dqDir, exist_ok=True)
regions = ['V1', 'ventral_stream', 'cortex']

# tSNR across all subjects (performed in anatomical space as native space registrations were poor)
print('calculating tSNR across all subjects...')
for region in regions:
    outDir = f'{dqDir}/tSNR/allSubjects/{region}'
    os.makedirs(outDir, exist_ok=True)
    plotData = {'topUp': {'b0': {'mean': 0,
                                 'std': 0},
                          'noB0': {'mean': 0,
                                 'std': 0}},
                'noTopUp': {'b0': {'mean': 0,
                                 'std': 0},
                          'noB0': {'mean': 0,
                                   'std': 0}}}
    for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
        for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
            means = []
            stds = []
            for subject in subjects:
                for session in experiment['scanInfo'][subject].keys():
                    sessDir = f'data/fMRI/individual/{subject}/{session}'
                    if 'restingState' in experiment['scanInfo'][subject][session]['funcScans'].keys():
                        dataPath = glob.glob(os.path.join(sessDir, f'functional/restingState/run01/inputs/rawData{topupString}{b0String}.nii*'))[0]
                        tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                        os.makedirs(tSNRdir, exist_ok=True)

                        # run motion correction
                        dataMotCor = os.path.join(tSNRdir, "motionCorrected.nii.gz")
                        if not os.path.exists(dataMotCor) or overwrite:
                            os.system(f'mcflirt -in {dataPath} -out {dataMotCor}')

                        # calculate tSNR
                        pathTmean = os.path.join(tSNRdir, "Tmean.nii.gz")
                        pathTstd = os.path.join(tSNRdir, "Tstd.nii.gz")
                        pathTSNR = os.path.join(tSNRdir, "tSNR.nii.gz")
                        if not os.path.exists(pathTmean) or overwrite:
                            print('Creating mean of timeseries...')
                            os.system(f'fslmaths {dataMotCor} -Tmean {pathTmean}')
                        if not os.path.exists(pathTstd) or overwrite:
                            print('Creating std of timeseries...')
                            os.system(f'fslmaths {dataMotCor} -Tstd {pathTstd}')
                        if not os.path.exists(pathTSNR) or overwrite:
                            print('Calculating tSNR map...')
                            os.system(f'fslmaths {pathTmean} -div {pathTstd} {pathTSNR}')

                        # convert masks from standard to func space
                        if region in ['V1', 'ventral_stream']:
                            funcMask = os.path.join(tSNRdir, f'{region}_func.nii.gz')
                            stdMask = glob.glob(f'/mnt/HDD12TB/masks/**/*{region}.nii.gz')[0]
                            std2funcMat = os.path.join(f'data/fMRI/individual/{subject}/{session}/functional/{regScan}/{regRun}/outputs/{topup}/{b0}/'
                                                    f'doubleGamma/firstLevel.feat/reg/standard2example_func.mat')  # get transformation matrix from a feat dir
                            exampleFunc = os.path.join(f'{os.getcwd()}/data/fMRI/individual/{subject}/{session}/functional/{regScan}/{regRun}/outputs/'
                                                       f'{topup}/{b0}/doubleGamma/firstLevel.feat/example_func.nii.gz')
                            if not os.path.exists(funcMask) or overwrite:
                                print('Converting standard space masks to native func space...')
                                os.system(f'flirt -in {stdMask} -ref {exampleFunc} -out {funcMask} -init {std2funcMat} -applyxfm -interp nearestneighbour')

                            # convert cortical mask to func space
                            ''' 
                            # method 1
                            regFile = f'{tSNRdir}/regFunc2Surf.lta'
                            if not os.path.isfile(regFile) or overwrite:
                                freesurferCommand = f'bbregister --bold --s {subject} --mov {pathTmean} --init-fsl --lta {regFile}'
                                os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                            bilatCortexMask = f'{tSNRdir}/cortex_bi.nii.gz'
                            if not os.path.exists(bilatCortexMask) or overwrite:
                                for hemi in ['lh', 'rh']:
                                    print(f'Converting {hemi} cortical mask...')
                                    inFile = f'/mnt/HDD12TB/freesurfer/subjects/{subject}/mri/{hemi}.ribbon.mgz'
                                    outFile = f'{tSNRdir}/cortex_{hemi}.nii.gz'
                                    freesurferCommand = f'mri_convert -at {regFile} -rt nearest {inFile} {outFile}'
                                    os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')
                                os.system(f'fslmaths {tSNRdir}/cortex_lh.nii.gz -add {tSNRdir}/cortex_rh.nii.gz {bilatCortexMask}')
                            # end of method 1 (does not give good results)
                            '''

                            # method 2
                            bilatCortexMask = f'{tSNRdir}/cortex_bi.nii.gz'
                            # convert freesurfer anat from mgz to nii
                            fsAnatFileMGZ = f'/mnt/HDD12TB/freesurfer/subjects/{subject}/mri/orig.mgz' # do not confuse with original nifti
                            fsAnatFileNII = f'{tSNRdir}/anatomical_fs.nii'
                            freesurferCommand = f'mri_convert --in_type mgz --out_type nii --out_orientation RAS {fsAnatFileMGZ} {fsAnatFileNII}'
                            os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                            # register freesurfer anat to orig anat (copy of anatomical saved by recon-all is not same as original nifti)
                            origAnatFile = f'/mnt/HDD12TB/freesurfer/subjects/{subject}/mri/anatomical.nii'
                            fsAnat2origAnatMat = f'{tSNRdir}/fsAnat2origAnat.mat'
                            os.system(f'flirt -in {fsAnatFileNII} -ref {origAnatFile} -omat {fsAnat2origAnatMat}')

                            # convert cortical mask to func space
                            for hemi in ['lh', 'rh']:
                                print(f'Converting {hemi} cortical mask...')

                                # cortex mgz to nifti
                                inFile = f'/mnt/HDD12TB/freesurfer/subjects/{subject}/mri/{hemi}.ribbon.mgz'
                                outFile = f'{tSNRdir}/cortex_{hemi}_FSanat.nii.gz'
                                if not os.path.exists(outFile):
                                    freesurferCommand = f'mri_convert --in_type mgz --out_type nii --out_orientation RAS {inFile} {outFile}'
                                    os.system(f'bash /mnt/HDD12TB/masterScripts/fMRI/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                # freesurfer anat space to orig anat space
                                inFile = outFile
                                outFile = f'{tSNRdir}/cortex_{hemi}_origAnat.nii.gz'
                                if not os.path.exists(outFile):
                                    os.system(f'flirt -in {inFile} -ref {origAnatFile} -applyxfm -init {fsAnat2origAnatMat} -interp nearestneighbour -out {outFile}')

                                # orig anat space to func space
                                inFile = outFile
                                outFile = f'{tSNRdir}/cortex_{hemi}.nii.gz'
                                origAnat2funcMat = os.path.join(f'data/fMRI/individual/{subject}/{session}/functional/{regScan}/{regRun}/outputs/{topup}/{b0}/'
                                                           f'doubleGamma/firstLevel.feat/reg/highres2example_func.mat')
                                os.system(f'flirt -in {inFile} -ref {pathTmean} -applyxfm -init {origAnat2funcMat} -interp nearestneighbour -out {outFile}')
                            os.system(f'fslmaths {tSNRdir}/cortex_lh.nii.gz -add {tSNRdir}/cortex_rh.nii.gz -bin {bilatCortexMask}')
                            # end of method 2

                            # combine cortex mask with ROI mask
                            finalMask = f'{tSNRdir}/{region}_cortex.nii.gz'
                            if not os.path.isfile(finalMask) or overwrite:
                                print('Combining ROI mask with cortical mask...')
                                os.system(f'fslmaths {tSNRdir}/cortex_bi.nii.gz -mul {funcMask} {finalMask}')

                        elif region == 'LGN':
                            finalMask = f'{tSNRdir}/LGN_bi.nii.gz'
                            leftHemi = os.path.join(sessDir, 'functional/figureGround_loc_v3/run03/outputs', topup, b0, 'doubleGamma/firstLevel.feat/masks/LGN_lh_00004_voxels.nii.gz')
                            rightHemi = os.path.join(sessDir, 'functional/figureGround_loc_v3/run03/outputs', topup, b0, 'doubleGamma/firstLevel.feat/masks/LGN_rh_00004_voxels.nii.gz')
                            os.system(f'fslmaths {leftHemi} -add {rightHemi} {finalMask}')

                        elif region == 'cortex':
                            finalMask = f'{tSNRdir}/cortex_bi.nii.gz'

                        means.append(float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -m').read()))
                        stds.append(float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -s').read()))

            plotData[topup][b0]['mean'] = np.mean(means)
            plotData[topup][b0]['std'] = np.mean(stds)

    b0_means = [plotData['topUp']['b0']['mean'], plotData['noTopUp']['b0']['mean']]
    b0_stds = [plotData['topUp']['b0']['std'], plotData['noTopUp']['b0']['std']]
    noB0_means = [plotData['topUp']['noB0']['mean'], plotData['noTopUp']['noB0']['mean']]
    noB0_stds = [plotData['topUp']['noB0']['std'], plotData['noTopUp']['noB0']['std']]

    x_pos = [0,1]
    fig, ax = plt.subplots()
    width = .33
    horzOffsets = [-.33, .33]
    fig.set_figwidth(6)
    ax.bar([-.167,1-.167], b0_means, width, yerr=b0_stds, align='center',
               ecolor='black', capsize=10, label='b0')
    ax.bar([.167,1.167], noB0_means, width, yerr=noB0_stds, align='center',
               ecolor='black', capsize=10, label='noB0')
    ax.set_ylabel('tSNR')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['topUp','noTopUp'])
    ax.set_title(f'{region} in all subjects')
    ax.yaxis.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{outDir}/tSNR_allSubjects_{region}.png')
    plt.show()
    plt.close()

# tSNR for each session of each subject
print('calculating tSNR for each subject and each session...')
for subject in subjects:
    for session in experiment['scanInfo'][subject].keys():
        sessDir = f'data/fMRI/individual/{subject}/{session}'
        if 'restingState' in experiment['scanInfo'][subject][session]['funcScans'].keys():
            for region in regions:
                outDir = f'{dqDir}/tSNR/{subject}/{session}/{region}'
                os.makedirs(outDir, exist_ok=True)
                plotData = {'topUp': {'b0': {'mean': 0,
                                             'std': 0},
                                      'noB0': {'mean': 0,
                                               'std': 0}},
                            'noTopUp': {'b0': {'mean': 0,
                                               'std': 0},
                                        'noB0': {'mean': 0,
                                                 'std': 0}}}
                for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
                    for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):

                        dataPath = glob.glob(os.path.join(sessDir, f'functional/restingState/run01/inputs/rawData{topupString}{b0String}.nii*'))[0]
                        tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                        pathTSNR = os.path.join(tSNRdir, 'tSNR.nii.gz')
                        finalMask = f'{tSNRdir}/{region}_cortex.nii.gz'
                        if region in ['cortex','LGN']:
                            finalMask = f'{tSNRdir}/{region}_bi.nii.gz'

                        # get tSNR values
                        plotData[topup][b0]['mean'] = float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -m').read())
                        plotData[topup][b0]['std'] = float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -s').read())

                b0_means = [plotData['topUp']['b0']['mean'], plotData['noTopUp']['b0']['mean']]
                b0_stds = [plotData['topUp']['b0']['std'], plotData['noTopUp']['b0']['std']]
                noB0_means = [plotData['topUp']['noB0']['mean'], plotData['noTopUp']['noB0']['mean']]
                noB0_stds = [plotData['topUp']['noB0']['std'], plotData['noTopUp']['noB0']['std']]

                x_pos = [0, 1]
                fig, ax = plt.subplots()
                width = .33
                horzOffsets = [-.33, .33]
                fig.set_figwidth(6)
                ax.bar([-.167, 1 - .167], b0_means, width, yerr=b0_stds, align='center',
                       ecolor='black', capsize=10, label='b0')
                ax.bar([.167, 1.167], noB0_means, width, yerr=noB0_stds, align='center',
                       ecolor='black', capsize=10, label='noB0')
                ax.set_ylabel('tSNR')
                ax.set_xticks(x_pos)
                ax.set_xticklabels(['topUp', 'noTopUp'])
                ax.set_title(f'resting state tSNR values in {subject} {session} {region}')
                ax.yaxis.grid(True)
                ax.legend()
                plt.tight_layout()
                plt.savefig(f'{outDir}/{subject}_{session}_{region}.png')
                plt.show()
                plt.close()

            # make histogram of tSNR values in cortical mask
            for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
                for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
                    tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                    pathTSNR = os.path.join(tSNRdir, 'tSNR.nii.gz')
                    corticalMask = f'{tSNRdir}/cortex_bi.nii.gz'
                    outDir = f'{dqDir}/tSNR/{subject}/{session}/cortex/{topup}/{b0}'
                    os.makedirs(outDir, exist_ok=True)
                    histData = [float(x) for x in os.popen(f'fslstats {pathTSNR} -k {corticalMask} -H 32 -8 248').read().split()]
                    max = float(os.popen(f'fslstats {pathTSNR} -k {corticalMask} -p 100').read())
                    min = float(os.popen(f'fslstats {pathTSNR} -k {corticalMask} -p 0').read())
                    x_pos = np.linspace(min,max,len(histData))
                    fig, ax = plt.subplots()
                    ax.bar(x_pos, histData, width=8)
                    ax.set_title(f'{subject}, {session}, {topup}, {b0} correction')
                    plt.xlabel('tSNR')
                    plt.ylabel('voxel count')
                    ax.yaxis.grid(True)
                    plt.xlim(-16,256)
                    plt.tight_layout()
                    plt.savefig(f'{outDir}/tSNR_hist_{topup}_{b0}.png')
                    plt.show()
                    plt.close()

            # make histogram of mean MR intensity values in cortical mask
            for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
                for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
                    tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                    pathTmean = os.path.join(tSNRdir, 'Tmean.nii.gz')
                    corticalMask = f'{tSNRdir}/cortex_bi.nii.gz'
                    outDir = f'{dqDir}/tSNR/{subject}/{session}/cortex/{topup}/{b0}'
                    os.makedirs(outDir, exist_ok=True)
                    max = float(os.popen(f'fslstats {pathTmean} -k {corticalMask} -p 100').read())
                    min = float(os.popen(f'fslstats {pathTmean} -k {corticalMask} -p 0').read())
                    histData = [float(x) for x in os.popen(f'fslstats {pathTmean} -k {corticalMask} -H 32 {min} {max}').read().split()]
                    x_pos = np.linspace(min, max, len(histData))
                    fig, ax = plt.subplots()
                    ax.bar(x_pos, histData, width = (max-min)/32)
                    ax.set_title(f'{subject}, {session}, {topup}, {b0} correction')
                    plt.xlabel('MR intensity')
                    plt.ylabel('voxel count')
                    ax.yaxis.grid(True)
                    plt.xlim(-10000, 400000)
                    plt.tight_layout()
                    plt.savefig(f'{outDir}/MRintensity_hist_{topup}_{b0}.png')
                    plt.show()
                    plt.close()

            # make scatterplot of tSNR values and signal intensity values in cortical mask
            for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
                for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):

                    tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                    outDir = f'{dqDir}/tSNR/{subject}/{session}/cortex/{topup}/{b0}'

                    pathTSNRanat = os.path.join(tSNRdir, 'tSNR.nii.gz')
                    pathTmeanAnat = f'{tSNRdir}/Tmean.nii.gz'
                    pathCorticalMask = f'{tSNRdir}/cortex_bi.nii.gz'

                    # apply mask
                    pathTSNRcortex = os.path.join(tSNRdir, 'tSNR_cortex.nii.gz')
                    pathTmeanCortex = f'{tSNRdir}/Tmean_cortex.nii.gz'
                    os.system(f'fslmaths {pathTSNRanat} -mul {pathCorticalMask} {pathTSNRcortex}')
                    os.system(f'fslmaths {pathTmeanAnat} -mul {pathCorticalMask} {pathTmeanCortex}')

                    # load masked data
                    dataTmean = nib.load(pathTmeanCortex).get_fdata().flatten()
                    dataTmean = np.ma.masked_where(dataTmean < MRthresh, dataTmean)
                    dataTSNR = nib.load(pathTSNRcortex).get_fdata().flatten()
                    dataTSNR = np.ma.masked_where(dataTmean < MRthresh, dataTSNR)
                    rVal = np.corrcoef(dataTSNR, dataTmean)[0,1]

                    # plot
                    fig, ax = plt.subplots()
                    ax.scatter(dataTmean, dataTSNR, s=.01)
                    ax.set_title(f'R = {rVal:.2f}, {topup}, {b0} correction')
                    plt.xlabel('signal intensity')
                    plt.ylabel('tSNR')
                    plt.xlim(-10000, 400000)
                    plt.ylim(-10, 200)
                    plt.tight_layout()
                    plt.savefig(f'{outDir}/tSNR_v_sigInt_{topup}_{b0}.png')
                    plt.show()
                    plt.close()

# make  plots for each top up / b0 combination
for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
    for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
        for region in regions:
            outDir = f'{dqDir}/tSNR/allSubjects/{region}/{topup}/{b0}'
            os.makedirs(outDir, exist_ok=True)
            means = []
            stds = []
            for subject in subjects:
                for session in experiment['scanInfo'][subject].keys():
                    sessDir = f'data/fMRI/individual/{subject}/{session}'
                    dataPath = glob.glob(os.path.join(sessDir, f'functional/restingState/run01/inputs/rawData{topupString}{b0String}.nii*'))[0]
                    tSNRdir = os.path.join(os.getcwd(), os.path.dirname(os.path.dirname(dataPath)), f'outputs/tSNR/{topup}/{b0}')
                    pathTSNR = os.path.join(tSNRdir, 'tSNR.nii.gz')
                    finalMask = f'{tSNRdir}/{region}_cortex.nii.gz'
                    if region in ['cortex','LGN']:
                        finalMask = f'{tSNRdir}/{region}_bi.nii.gz'

                    # get tSNR values
                    means.append(float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -m').read()))
                    stds.append(float(os.popen(f'fslstats {pathTSNR} -k {finalMask} -s').read()))

            x_pos = np.arange(nSubjects)
            fig, ax = plt.subplots()
            width = .8
            fig.set_figwidth(6)
            ax.bar(x_pos, means, width, yerr=stds, align='center',
                   ecolor='black', capsize=10, label='b0')
            ax.set_ylabel('tSNR')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(subjects)
            ax.set_title(f'resting state tSNR values in {region} {topup} {b0}')
            ax.yaxis.grid(True)
            ax.legend()
            plt.tight_layout()
            plt.savefig(f'{outDir}/{region}_{topup}_{b0}.png')
            plt.show()
            plt.close()