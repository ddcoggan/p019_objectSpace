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
cropParams = {'F013': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'F016': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'F111': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'F118': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'F120': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'F121': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'M012': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'M015': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}},
              'M119': {'lh': {'coords': [0,175],
                              'size': [250,415]},
                       'rh': {'coords': [0,175],
                              'size': [250,415]}}}

subjects = list(cropParams.keys())

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

                    contrastDirs = sorted(glob.glob(os.path.join('analysis/results/fMRI/surfacePlots', topup,
                                                                 b0, HRFmodel, subject, scan, '*')))[2:]
                    contrastDirs += sorted(glob.glob(os.path.join('analysis/results/fMRI/surfacePlots', topup,
                                                                 b0, HRFmodel, subject, scan, 'all/*/non_overla')))

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

