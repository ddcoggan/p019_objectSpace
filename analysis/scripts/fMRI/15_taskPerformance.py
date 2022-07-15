
import sys
import os
import glob
import pandas as pd
import pickle
import itertools
import numpy as np
from scipy.io import loadmat
from argparse import Namespace
import matplotlib.pyplot as plt
import datetime
import shutil

# get scan info from experiment file
from analysis.scripts.fMRI.experiment import experiment

taskData = {'subject': [],
            'session': [],
            'scan': [],
            'run': [],
            'accuracy': [],
            'RT': []}
for subject in experiment['scanInfo'].keys():
    for s, session in enumerate(experiment['scanInfo'][subject].keys()):
        print(f'Subject: {subject}, Session: {session}')

        sessID = experiment['scanInfo'][subject][session]['sessID']
        sessDir = os.path.join('data/fMRI/individual', subject, session)

        for scan in experiment['scanInfo'][subject][session]['funcScans'].keys():
            for r, run in enumerate(experiment['scanInfo'][subject][session]['funcScans'][scan]):

                # generate event files from log files
                if scan != 'restingState':
                    logFile = glob.glob(os.path.join('data/fMRI/events/logFiles', subject, session, scan, f'*_rn{r+1}_*'))[0]
                    logData = loadmat(logFile)
                    responses = logData['experiment']['targetResponseTimes'].item()
                    if len(responses) > 0:
                        accuracy = len(responses[0])/32 # corrects for error in stim pres code which assumed 40 trials (not 32)
                        RT = 0.667-np.mean(responses[0])
                    else:
                        accuracy = 0
                        RT = None
                    taskData['subject'].append(subject)
                    taskData['session'].append(session)
                    taskData['scan'].append(scan)
                    taskData['run'].append(r+1)
                    taskData['accuracy'].append(accuracy)
                    taskData['RT'].append(RT)

taskDF = pd.DataFrame(taskData)
