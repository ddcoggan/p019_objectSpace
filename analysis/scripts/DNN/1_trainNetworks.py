'''
uses train.py from masterScripts/DNN to train different DNNs on different datasets
if script is interrupted or additional epochs are required at a later time, script will continue from last recorded epoch
'''

import os
import glob
import sys

sys.path.append('/mnt/HDD12TB/masterScripts/DNN')
from train import train

overwrite = False
models = ['alexnet']
datasets = ['imagenet1000_downsampled']#, 'places365_standard']
learningRate = .01
optimizer = 'SGD'
batchSize = 2048
nEpochs = 65
workers = 8

for model in models:

    for dataset in datasets:

        outDir = os.path.join(os.getcwd(), 'data', model, dataset)

        # get restart from file if necessary
        weightFiles = sorted(glob.glob(os.path.join(outDir, 'params/*.pt')))
        if 0 < len(weightFiles):
            restartFrom = weightFiles[-1]
        else:
            restartFrom = None

        # call script
        if len(weightFiles) < nEpochs+1 or overwrite:
            train(model, dataset, False, learningRate, optimizer, batchSize, nEpochs, restartFrom, workers, outDir, None)
