import os
import sys
import glob
import datetime
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

regionOrder = ['eba','fba','ofa','ffa','psts','opa','ppa','rsc']
NMLregionOrder = ['nml_lat_1','nml_lat_2','nml_med_1','nml_med_2']
catNames = ['body','face','house','object']
condNames = ['animate_spiky', 'animate_stubby', 'inanimate_spiky', 'inanimate_stubby']

dataFile = 'data/fMRI/group/overlapROIknownRegions.csv'
data = pd.read_csv(dataFile)
outDir = 'analysis/results/fMRI/overlapROI/knownRegions'
os.makedirs(outDir, exist_ok=True)
colours = ['chartreuse', 'deepskyblue', 'orange', 'fuchsia']
'''
for size in ['16','64','256','allVoxels']:
    for c, cond in enumerate(condNames):
        theseData = data[(data['nVoxTarget'] == size) & (data['condition'] == cond)]
        df = theseData.groupby(['region','voxType']).agg(['mean']).reset_index().iloc[:, [0, 1, 3]]
        df.columns = df.columns.droplevel(1)
        df.index = df['region']
        df = df.pivot(index='region', columns='voxType', values='nVox').reindex(index = regionOrder)
        dfPlot = df[['regionUnique','overlap','conditionUnique']].plot.bar(stacked=True, color=['gray','black',colours[c]])
        fig = dfPlot.get_figure()
        plt.ylabel('voxels')
        plt.xticks(rotation = 'horizontal')
        plt.legend(['ROI','overlap', f'{condNames[c]} cluster'], loc='center left', bbox_to_anchor=(1.0, 0.5)).set_title('')
        plt.tight_layout()
        fig.savefig(os.path.join(outDir, f'{cond}_{size}.png'))
        plt.show(block=False)

'''
dataFile = 'data/fMRI/group/overlapROINMLRegions.csv'
data = pd.read_csv(dataFile)
outDir = 'analysis/results/fMRI/overlapROI/NMLRegions'
os.makedirs(outDir, exist_ok=True)
colours = ['grey','red','deepskyblue','chartreuse']
conds = ['unique','house','face','body']
legendLabels = ['body-selective','face-selective','scene-selective','no other selectivity']
legendLabels.reverse()
for size in ['16','64','256','allVoxels']:
    theseData = data[(data['nVoxTarget'] == size)]
    if size == 'allVoxels':
        df = theseData.groupby(['region', 'condition']).agg(['mean']).reset_index().iloc[:, [0, 1, 3]]
    else:
        df = theseData.groupby(['region', 'condition']).agg(['mean']).reset_index().iloc[:, [0, 1, 4]]

    df.columns = df.columns.droplevel(1)
    df.index = df['region']
    df = df.pivot(index='region', columns='condition', values='nVox').reindex(index = NMLregionOrder)
    dfPlot = df[conds].plot.bar(stacked=True, color=colours,figsize=(7*0.7,4*0.7))
    fig = dfPlot.get_figure()
    plt.ylabel('number of voxels')
    plt.xlabel(None)
    plt.xticks(np.arange(4), rotation = 'horizontal', labels=['NML\nlat 1','NML\nlat 2','NML\nmed 1','NML\nmed 2'])
    plt.legend(labels = legendLabels,bbox_to_anchor=(1.0, 0.7), loc='upper left', labelspacing = -2.5, frameon=False).set_title('')
    #plt.tight_layout()
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0+.05, box.y0+.1, box.width * 0.55, box.height])
    fig.savefig(os.path.join(outDir, f'{size}.png'))
    plt.show(block=False)
'''
dataFile = 'data/fMRI/group/overlapROIallClusters.csv'
data = pd.read_csv(dataFile)
outDir = 'analysis/results/fMRI/overlapROI/allClusters'
os.makedirs(outDir, exist_ok=True)
for c, cond in enumerate(condNames):
    theseData = data[data['AARcond'] == cond]
    df = theseData.groupby(['BFHOcond', 'voxType']).agg(['mean']).reset_index().iloc[:, [0, 1, 3]]
    df.columns = df.columns.droplevel(1)
    df.index = df['BFHOcond']
    df = df.pivot(index='BFHOcond', columns='voxType', values='nVox')
    dfPlot = df[['BFHOunique', 'overlap', 'AARunique']].plot.bar(stacked=True,color=['gray', 'black', colours[c]])
    fig = dfPlot.get_figure()
    plt.ylabel('voxels')
    plt.xticks(rotation='horizontal')
    plt.legend(['category cluster', 'overlap', f'{condNames[c]} cluster'], loc='center left', bbox_to_anchor=(1.0, 0.5)).set_title('')
    plt.tight_layout()
    fig.savefig(os.path.join(outDir, f'{cond}.png'))
    plt.show(block=False)
'''