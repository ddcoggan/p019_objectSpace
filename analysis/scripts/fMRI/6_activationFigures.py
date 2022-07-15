import os
import glob
import datetime
from PIL import Image
from PIL import ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True
import numpy as np
from analysis.scripts.fMRI.experiment import experiment

overwrite = 1
hemis = ['lh', 'rh']
views = ['inferior', 'lateral']
imSize = {'inferior': [235, 415],
          'lateral': [300, 415]}
cropParams = {'inferior': {'lh': [0, 175],
                           'rh': [0, 175]},
              'lateral': {'lh': [290, 0],
                          'rh': [0, 0]}}
subjects = list(experiment['scanInfo'].keys())
'''
#subjects.append('allSubjects')
for subject in subjects:
    if subject != 'allSubjects':
        session = list(experiment['scanInfo'][subject].keys())[0]
    else:
        session = 'allSessions'
    for scan in experiment['design'].keys():
        for topup in ['topUp']: #, 'noTopUp']:
            for b0 in ['noB0']: # 'b0']:
                for HRFmodel in ['doubleGamma']:#, 'singleGamma']:

                    print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | topup: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                    contrastDirs = sorted(glob.glob(os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, '????*')))

                    for contrastDir in contrastDirs:
                        for mapType in ['cluster','zstat']:

                            finalFile = os.path.join(contrastDir, mapType, 'allViews.png')
                            testFile = os.path.join(contrastDir, mapType, 'inferior_rh.png')
                            if os.path.isfile(testFile):
                                if not os.path.isfile(finalFile) or overwrite:

                                    # load and crop images
                                    images = {'inferior': {'lh': [],
                                                          'rh': []},
                                              'lateral': {'lh': [],
                                                          'rh': []}}
                                    for view in views:
                                        for hemi in hemis:
                                            imPath = os.path.join(contrastDir, mapType, f'{view}_{hemi}.png')
                                            imOrig = Image.open(imPath)
                                            cropBox = np.array(cropParams[view][hemi])
                                            cropBox = np.append(cropBox, cropBox+imSize[view])
                                            imCropped = imOrig.crop(cropBox)
                                            images[view][hemi] = imCropped
                                    concatIm = Image.new('RGB', (imSize['lateral'][0]*2 + imSize['inferior'][0]*2, imSize['lateral'][1]))
                                    concatIm.paste(images['lateral']['rh'], (0, 0))
                                    concatIm.paste(images['inferior']['rh'], (imSize['lateral'][0], 0))
                                    concatIm.paste(images['inferior']['lh'], (imSize['lateral'][0] + imSize['inferior'][0], 0))
                                    concatIm.paste(images['lateral']['lh'], (imSize['lateral'][0] + imSize['inferior'][0]*2, 0))
                                    concatIm.save(finalFile)

                        for contrast in ['versusOthers','noContrast']:
                            for mapType in ['cluster', 'zstat']:
                                contrastDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', mapType, contrast)

                                finalFile = os.path.join(contrastDir, 'allViews.png')
                                testFile = os.path.join(contrastDir, 'inferior_rh.png')
                                if os.path.isfile(testFile):
                                    if not os.path.isfile(finalFile) or overwrite:

                                        # load and crop images
                                        images = {'inferior': {'lh': [],
                                                              'rh': []},
                                                  'lateral': {'lh': [],
                                                              'rh': []}}
                                        for view in views:
                                            for hemi in hemis:
                                                imPath = os.path.join(contrastDir, f'{view}_{hemi}.png')
                                                imOrig = Image.open(imPath)
                                                cropBox = np.array(cropParams[view][hemi])
                                                cropBox = np.append(cropBox, cropBox + imSize[view])
                                                imCropped = imOrig.crop(cropBox)
                                                images[view][hemi] = imCropped
                                        concatIm = Image.new('RGB', (
                                        imSize['lateral'][0] * 2 + imSize['inferior'][0] * 2, imSize['lateral'][1]))
                                        concatIm.paste(images['lateral']['rh'], (0, 0))
                                        concatIm.paste(images['inferior']['rh'], (imSize['lateral'][0], 0))
                                        concatIm.paste(images['inferior']['lh'], (imSize['lateral'][0] + imSize['inferior'][0], 0))
                                        concatIm.paste(images['lateral']['lh'],
                                                       (imSize['lateral'][0] + imSize['inferior'][0] * 2, 0))
                                        concatIm.save(finalFile)
'''
# repeat for just inferior views
imSize = [230, 375]
cropParams = {'lh': [2, 200],
              'rh': [0, 200]}

#subjects.append('allSubjects')
for subject in subjects:
    if subject != 'allSubjects':
        session = list(experiment['scanInfo'][subject].keys())[0]
    else:
        session = 'allSessions'
    for scan in experiment['design'].keys():
        for topup in ['topUp']: #, 'noTopUp']:
            for b0 in ['noB0']: # 'b0']:
                for HRFmodel in ['doubleGamma']:#, 'singleGamma']:

                    print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | Session: {session} | Scan: {scan} | topup: {topup} | b0: {b0} | HRF model: {HRFmodel}')

                    '''
                    contrastDirs = sorted(glob.glob(os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, '????*')))

                    for contrastDir in contrastDirs:
                        for mapType in ['cluster','zstat']:

                            finalFile = os.path.join(contrastDir, mapType, 'inferiorViews.png')
                            testFile = os.path.join(contrastDir, mapType, 'inferior_rh.png')
                            if os.path.isfile(testFile):
                                if not os.path.isfile(finalFile) or overwrite:


                                    # load and crop images
                                    images = {'lh': [],
                                              'rh': []}
                                    for hemi in hemis:
                                        imPath = os.path.join(contrastDir, mapType, f'inferior_{hemi}.png')
                                        imOrig = Image.open(imPath)
                                        cropBox = np.array(cropParams[hemi])
                                        cropBox = np.append(cropBox, cropBox+imSize)
                                        imCropped = imOrig.crop(cropBox)
                                        images[hemi] = imCropped
                                    concatIm = Image.new('RGB', (imSize[0]*2, imSize[1]))
                                    concatIm.paste(images['rh'], (0, 0))
                                    concatIm.paste(images['lh'], (imSize[0], 0))
                                    concatIm.save(finalFile)

                        for contrast in ['versusOthers','noContrast']:
                            for mapType in ['cluster', 'zstat']:
                                contrastDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan, 'all', mapType, contrast)

                                finalFile = os.path.join(contrastDir, 'inferiorViews.png')
                                testFile = os.path.join(contrastDir, 'inferior_rh.png')
                                if os.path.isfile(testFile):
                                    if not os.path.isfile(finalFile) or overwrite:

                                        # load and crop images
                                        images = {'lh': [],
                                                  'rh': []}
                                        for hemi in hemis:
                                            imPath = os.path.join(contrastDir, f'inferior_{hemi}.png')
                                            imOrig = Image.open(imPath)
                                            cropBox = np.array(cropParams[hemi])
                                            cropBox = np.append(cropBox, cropBox + imSize)
                                            imCropped = imOrig.crop(cropBox)
                                            images[hemi] = imCropped
                                        concatIm = Image.new('RGB', (imSize[0] * 2, imSize[1]))
                                        concatIm.paste(images['rh'], (0, 0))
                                        concatIm.paste(images['lh'], (imSize[0], 0))
                                        concatIm.save(finalFile)
                    '''
# PPA inanimate spiky > inanimate stubby
imSize = [210, 500]
cropParams = {'lh': [480, 60],
              'rh': [295, 60]}
overwrite = 1

# subjects.append('allSubjects')
for subject in subjects:
    if subject != 'allSubjects':
        session = list(experiment['scanInfo'][subject].keys())[0]
    else:
        session = 'allSessions'
    for scan in experiment['design'].keys():
        for topup in ['topUp']:  # , 'noTopUp']:
            for b0 in ['noB0']:  # 'b0']:
                for HRFmodel in ['doubleGamma']:  # , 'singleGamma']:

                    if scan == 'animacyAspectRatioLoc_v2':
                        contrastDir = os.path.join('analysis/results/fMRI/surfacePlots', topup, b0, HRFmodel, subject, scan,
                                         'inanimate_spiky>inanimate_stubby')

                        finalFile = os.path.join(contrastDir, 'copeInPPA_InferiorViews.png')
                        testFile = os.path.join(contrastDir, 'copeInPPA_rh_inferior.png')
                        if os.path.isfile(testFile):
                            if not os.path.isfile(finalFile) or overwrite:

                                # load and crop images
                                images = {'lh': [],
                                          'rh': []}
                                for hemi in hemis:
                                    imPath = os.path.join(contrastDir,f'copeInPPA_{hemi}_inferior.png')
                                    imOrig = Image.open(imPath)
                                    cropBox = np.array(cropParams[hemi])
                                    cropBox = np.append(cropBox, cropBox + imSize)
                                    imCropped = imOrig.crop(cropBox)
                                    images[hemi] = imCropped
                                concatIm = Image.new('RGB', (imSize[0] * 2, imSize[1]))
                                concatIm.paste(images['rh'], (0, 0))
                                concatIm.paste(images['lh'], (imSize[0], 0))
                                concatIm.save(finalFile)

