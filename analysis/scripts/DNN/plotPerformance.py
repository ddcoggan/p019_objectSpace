'''
makes plots similar to wandb for trained networks.
Note: basic plots are automatically created during training but can be regenerated
using log files with desired style using this script
'''

import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import pickle

netDirs = sorted(glob.glob('data/*/*'))

for netDir in netDirs:
    
    print(netDir)
    plotDir = os.path.join(netDir, 'plots')
    os.makedirs(plotDir, exist_ok=True)
    
    data = pd.read_csv(os.path.join(netDir, 'log.csv'))
    
    nEpochs = max(data['epoch'])
    epochs = list(range(nEpochs+1))
    

    for epoch in epochs:
        if epoch == 0:
            theseData = data[data['epoch'] == epoch]
            performance = {'train': {'acc1': [None], 'acc5': [None], 'loss': [None]},
                           'eval': {'acc1': [theseData.iloc[-1]['cumAcc1epoch']], 'acc5': [theseData.iloc[-1]['cumAcc5epoch']], 'loss': [theseData.iloc[-1]['cumLossEpoch']]}}
        else:
            for trainEval in ['train', 'eval']:

                theseData = data[(data['epoch'] == epoch) & (data['trainEval'] == trainEval)]

                performance[trainEval]['acc1'].append(theseData.iloc[-1]['cumAcc1epoch'])
                performance[trainEval]['acc5'].append(theseData.iloc[-1]['cumAcc5epoch'])
                performance[trainEval]['loss'].append(theseData.iloc[-1]['cumLossEpoch'])
    
    for plotType in ['acc1', 'acc5', 'loss']:
        
        plt.plot(epochs, performance['train'][plotType], label='train')
        plt.plot(epochs, performance['eval'][plotType], label='test')
        plt.xlabel('epoch')
        plt.ylabel(plotType)
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(plotDir, f'{plotType}.png'))
        plt.show()
        plt.close()

    plotStatsFile = os.path.join(netDir, 'plotStats.pkl')
    pickle.dump(performance, open(plotStatsFile, 'wb'))
        