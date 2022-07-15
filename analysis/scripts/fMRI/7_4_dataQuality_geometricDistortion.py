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

newSubjects = 1 # overwrites outputted data file if 1

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

def calc_MI(X,Y,bins):

    c_XY = np.histogram2d(X,Y,bins)[0] # 2d array of counts
    c_X = np.histogram(X,bins)[0] # 1d array of counts
    c_Y = np.histogram(Y,bins)[0] # 1d array of counts

    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)

    MI = (H_X + H_Y) / H_XY
    return MI

def shan_entropy(c):
    c_normalized = c / float(np.sum(c)) # represent count as proportion of total
    c_normalized = c_normalized[np.nonzero(c_normalized)] # remove zeros
    H = -sum(c_normalized * np.log2(c_normalized)) # get log2 transform, multiply with original, sum values, invert
    return H

# store collated data in a dictionary of lists
maxSessions = 1
subjects = list(experiment['scanInfo'].keys())
nSubjects = len(subjects)
dqDir = 'analysis/results/fMRI/dataQuality'
os.makedirs(dqDir, exist_ok=True)
dataFile = f'{dqDir}/geometricDistortion/mutualInfo_R2.pkl'

# normalised mutual information and R2
os.makedirs(f'{dqDir}/geometricDistortion', exist_ok=True)
if newSubjects:
    data = {}
    tempDir = f'{dqDir}/geometricDistortion/temp'  # set up temp dir
    os.makedirs(tempDir, exist_ok=True)
    for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
        data[topup] = {}
        for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
            data[topup][b0] = {}
            for subject in experiment['scanInfo'].keys():
                data[topup][b0][subject] = {'betweenSession': {'mutualInfo': [],
                                                               'R2': []}}
                imageData = {} # store images across sessions
                nSessions = len(list(experiment['scanInfo'][subject].keys()))
                maxSessions = max(maxSessions, nSessions)
                for session in experiment['scanInfo'][subject].keys():
                    sessDir = f'data/fMRI/individual/{subject}/{session}'
                    data[topup][b0][subject][session] = {'betweenVol': {'mutualInfo': [], 'R2': []},
                                                         'betweenRun': {'mutualInfo': [], 'R2': []},
                                                         'versusAnatomical': {'mutualInfo': [], 'R2': []},
                                                         'versusFuncNoEPI': {'mutualInfo': [], 'R2': []}}

                    imageData[session] = []  # list for first volume to compare between sessions

                    for scan in experiment['design'].keys():
                        for run in range(len(experiment['scanInfo'][subject][session]['funcScans'][scan])):

                            '''
                            # for debugging
                            topup = 'topUp'
                            topupString = '_topUp'
                            b0 = 'b0'
                            b0String = '_b0'
                            subject = 'M015'
                            session = '201215'
                            scan = 'figureGround_v3'
                            run = 0
                            '''

                            print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | topUp: {topup} | b0: {b0} | subject: {subject} | session: {session} | scan: {scan} | run: {run+1}')

                            # BETWEEN VOLUME
                            comparison = 'betweenVolume'
                            dataPath = glob.glob(os.path.join(sessDir, f'functional/{scan}/run{run+1:02}/inputs/rawData{topupString}{b0String}.nii*'))[0]
                            firstVol = nib.load(dataPath).get_fdata()[:, :, :, 0].flatten()
                            lastVol = nib.load(dataPath).get_fdata()[:, :, :, -1].flatten()
                            data[topup][b0][subject][session]['betweenVol']['mutualInfo'].append(calc_MI(firstVol, lastVol, bins=64))
                            data[topup][b0][subject][session]['betweenVol']['R2'].append(np.corrcoef(firstVol, lastVol)[0, 1] ** 2)

                            # store first volume for between run and session analysis
                            imageData[session].append(firstVol)

                            # VERSUS ANATOMICAL
                            anatomicalPath = os.path.join(sessDir, 'anatomical/anatomical.nii')  # get path to anatomical
                            anatCopyPath = os.path.join(tempDir, 'anatomical.nii')  # new path for anatomical in temp dir
                            shutil.copy(anatomicalPath, anatCopyPath)  # copy anatomical over
                            dataTmeanPath = os.path.join(tempDir, 'funcTmean.nii.gz')  # set path for average func image
                            os.system(f'fslmaths {dataPath} -Tmean {dataTmeanPath}')  # calculate average func image
                            transMat = os.path.join(tempDir, 'anat_to_func.mat')  # make transformation matrix from highres to func space
                            os.system(f'flirt -in {anatCopyPath} -ref {dataTmeanPath} -omat {transMat}')  # make transmat if necessary
                            funcAnat = os.path.join(tempDir,'anatInFuncSpace.nii.gz')  # path to anatomical volume in func space
                            os.system(f'flirt -in {anatCopyPath} -ref {dataTmeanPath} -out {funcAnat} -init {transMat} -applyxfm -interp trilinear')  # transform anatomical
                            funcMask = os.path.join(tempDir, 'funcMask.nii.gz')
                            os.system(f'fslmaths {dataTmeanPath} -bin {funcMask}')  # make binary mask of func data
                            funcAnatMasked = os.path.join(tempDir, 'anatInFuncSpaceMasked.nii.gz')
                            os.system(f'fslmaths {funcAnat} -mul {funcMask} {funcAnatMasked}')  # apply mask to anatomical in func space
                            anatVol = -nib.load(funcAnatMasked).get_fdata().flatten()  # load anatomical volume and invert
                            data[topup][b0][subject][session]['versusAnatomical']['mutualInfo'].append(calc_MI(firstVol, anatVol, bins=64))
                            data[topup][b0][subject][session]['versusAnatomical']['R2'].append(np.corrcoef(firstVol, anatVol)[0, 1] ** 2)


                            # VERSUS FUNC NO EPI
                            funcNoEPIpath = os.path.join(sessDir, 'funcNoEPI/funcNoEPI.nii')  # get path to anatomical
                            funcNoEPIvol = -nib.load(funcNoEPIpath).get_fdata().flatten()  # load volume and invert
                            data[topup][b0][subject][session]['versusFuncNoEPI']['mutualInfo'].append(calc_MI(firstVol, funcNoEPIvol, bins=64))
                            data[topup][b0][subject][session]['versusFuncNoEPI']['R2'].append(np.corrcoef(firstVol, funcNoEPIvol)[0, 1] ** 2)

                    # BETWEEN RUN (WITHIN SESSION)
                    for a in range(len(imageData[session])):
                        for b in range(a + 1, len(imageData[session])):
                            data[topup][b0][subject][session]['betweenRun']['mutualInfo'].append(calc_MI(imageData[session][a], imageData[session][b], bins=64))
                            data[topup][b0][subject][session]['betweenRun']['R2'].append(np.corrcoef(imageData[session][a], imageData[session][b])[0, 1] ** 2)

                # BETWEEN SESSION
                sessions = list(experiment['scanInfo'][subject].keys())
                if len(sessions) == 2:
                    for a in range(len(imageData[sessions[0]])):
                        for b in range(len(imageData[sessions[1]])):
                            data[topup][b0][subject]['betweenSession']['mutualInfo'].append(calc_MI(imageData[sessions[0]][a], imageData[sessions[1]][b], bins=64))
                            data[topup][b0][subject]['betweenSession']['R2'].append(np.corrcoef(imageData[sessions[0]][a], imageData[sessions[1]][b])[0, 1] ** 2)

    pkl.dump(data, open(dataFile, 'wb'))
    shutil.rmtree(tempDir)
else:
    data = pkl.load(open(dataFile, 'rb'))

# make plots
ybounds = {'mutualInfo': [0,2], 'R2': [0,1]}
figwidth = .8
for topup in ['topUp', 'noTopUp']:
    for b0 in ['b0', 'noB0']:

        outDir = os.path.join(dqDir, 'geometricDistortion', topup, b0)
        os.makedirs(outDir, exist_ok=True)

        # BETWEEN SESSION
        for measure, measureString in zip(['mutualInfo','R2'],['normalised mutual information', 'R2']):
            plotData = [[],[]]
            MSsubs = []
            for subject in experiment['scanInfo'].keys():
                sessions = list(experiment['scanInfo'][subject].keys())
                if len(sessions) == 2:
                    MSsubs.append(subject)
                    plotData[0].append(np.mean(data[topup][b0][subject]['betweenSession'][measure]))
                    plotData[1].append(sem(data[topup][b0][subject]['betweenSession'][measure]))

            nMSsubjects = len(MSsubs) # number of subject with multiple sessions
            x_pos = np.arange(nMSsubjects)
            fig, ax = plt.subplots()
            fig.set_figwidth(figwidth+nMSsubjects*.5)
            ax.bar(x_pos, plotData[0], yerr=plotData[1], align='center', ecolor='black', capsize=10)
            ax.set_ylabel(measure)
            ax.set_xlabel('subject')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(MSsubs)
            ax.set_title(f'Between session {measureString}, {topup}, {b0}')
            ax.yaxis.grid(True)
            ax.set_ybound(ybounds[measure])
            plt.tight_layout()
            plt.savefig(f'{outDir}/betweenSession_{measure}.png')
            plt.show()
            plt.close()

        # BETWEEN RUN
        for measure, measureString in zip(['mutualInfo', 'R2'], ['normalised mutual information', 'R2']):
            plotData = []
            for sess in range(maxSessions):
                toAppend = [[],[]]
                for subject in experiment['scanInfo'].keys():
                    sessions = list(experiment['scanInfo'][subject].keys())
                    if sess < len(sessions):
                        toAppend[0].append(np.mean(data[topup][b0][subject][sessions[sess]]['betweenRun'][measure]))
                        toAppend[1].append(sem(data[topup][b0][subject][sessions[sess]]['betweenRun'][measure]))
                    else:
                        toAppend[0].append(0)
                        toAppend[1].append(0)
                plotData.append(toAppend)

            x_pos = np.arange(nSubjects)
            fig, ax = plt.subplots()
            width = 1/(maxSessions+1)
            horzOffsets = np.linspace(width, 1-width, maxSessions)-.5
            fig.set_figwidth(figwidth+nSubjects*.5)
            for sess in range(maxSessions):
                ax.bar(x_pos + horzOffsets[sess], plotData[sess][0], width, yerr=plotData[sess][1], align='center', ecolor='black', capsize=10, label = f'session {sess+1}')
            ax.set_ylabel(measure)
            ax.set_xlabel('subject')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(subjects)
            ax.set_title(f'Between run {measureString}, {topup}, {b0}')
            ax.yaxis.grid(True)
            ax.set_ybound(ybounds[measure])
            #ax.legend()
            plt.tight_layout()
            plt.savefig(f'{outDir}/betweenRun_{measure}.png')
            plt.show()
            plt.close()

        # WITHIN RUN COMPARISONS
        for comparison in ['betweenVol','versusAnatomical','versusFuncNoEPI']:
            for measure, measureString in zip(['mutualInfo', 'R2'], ['normalised mutual information', 'R2']):
                plotData = []
                for sess in range(maxSessions):
                    toAppend = [[], []]
                    for subject in experiment['scanInfo'].keys():
                        sessions = list(experiment['scanInfo'][subject].keys())
                        if sess < len(sessions):
                            toAppend[0].append(np.mean(data[topup][b0][subject][sessions[sess]][comparison][measure]))
                            toAppend[1].append(sem(data[topup][b0][subject][sessions[sess]][comparison][measure]))
                        else:
                            toAppend[0].append(0)
                            toAppend[1].append(0)
                    plotData.append(toAppend)

                x_pos = np.arange(nSubjects)
                fig, ax = plt.subplots()
                width = 1 / (maxSessions + 1)
                horzOffsets = np.linspace(width, 1 - width, maxSessions) - .5
                fig.set_figwidth(figwidth+nSubjects*.5)
                for sess in range(maxSessions):
                    ax.bar(x_pos + horzOffsets[sess], plotData[sess][0], width, yerr=plotData[sess][1], align='center', ecolor='black', capsize=10, label = f'session {sess+1}')
                ax.set_ylabel(measure)
                ax.set_xlabel('subject')
                ax.set_xticks(x_pos)
                ax.set_xticklabels(subjects)
                ax.set_title(f'{comparison} {measureString}, {topup}, {b0}')
                ax.yaxis.grid(True)
                #ax.legend()
                ax.set_ybound(ybounds[measure])
                plt.tight_layout()
                plt.savefig(f'{outDir}/{comparison}_{measure}.png')
                plt.show()
                plt.close()
