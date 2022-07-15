"""
plots PCA space in neural networks
"""
import os
import sys
import glob
import datetime
import pickle
import matplotlib.pyplot as plt
import numpy as np
import torch
import shutil
from sklearn.decomposition import PCA
from PIL import Image
from PIL import ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True
overwrite = 1
sys.path.append('/mnt/HDD12TB/masterScripts/DNN')
from saveOutputs import saveOutputs
from modelLayerLabels import modelLayerLabels
from centreCropResize import centreCropResize

### CONFIGURATION
models = ['alexnet']
databases = {'train': ['imagenet1000'],#'places365_standard'],
             'components': ['BaoTsao', 'BaoTsao2', 'imagenet1000_subset', 'imagenet1000_subset2'],#, 'BaoTsao2', 'imagenet1000_subset', 'places365_subset', 'BOSS', 'BigSmallObjects'],
             'samples': ['ecoset_subset', 'BaoTsao3', 'BFHO', 'imagenet1000_subset', 'imagenet1000_subset2']}#, 'imagenet1000_subset', 'places365_subset', 'BOSS', 'BigSmallObjects', 'BaoTsao', 'BaoTsao2']}
layers2use = ['fc6'] # any layer starting with these strings is used
weights = 'final' # select which weights to use. 'final' is last training epoch. TODO 'maxEval' is epoch with highest evaluation performance.
nPCs = 8
nPCAimages = 16

for model in models:

    layers = modelLayerLabels[model]
    theseLayers = []
    for layer in layers:
        for layer2use in layers2use:
            if layer.startswith(layer2use):
                theseLayers.append(layer)

    for dbt in databases['train']:

        modelDir = os.path.join('data/DNN', model, dbt)
        paramsDir = os.path.join(modelDir, 'params')
        responseDir = os.path.join(modelDir, 'responses')

        if weights == 'final':
            #paramsFile = sorted(glob.glob(os.path.join(paramsDir, '*.pt')))[-1]
            paramsFile = None

        # measure responses
        for dbr in databases['components'] + databases['samples']:

            # measure responses
            thisResponseDir = os.path.join(responseDir, dbr)
            os.makedirs(thisResponseDir, exist_ok=True)
            imageDir = os.path.join('/home/dave/Datasets', dbr)
            if not os.path.isdir(imageDir):
                imageDir = os.path.join('/media/dave/HotProjects/Datasets', dbr)
            imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.png'), recursive=True))
            if len(imageFiles) == 0:
                imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.jpg'), recursive=True))
            if len(imageFiles) == 0:
                imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.tif'), recursive=True))
            if len(imageFiles) == 0:
                raise Exception('No png, jpg or tif image files found!')

            # loop over images and measure responses
            doneFiles = []
            for i, imageFile in enumerate(imageFiles):
                outFile = f'{thisResponseDir}/{i + 1:08}_{os.path.basename(imageFile).split(sep=".")[0]}.pkl'
                if not os.path.isfile(outFile) and imageFile:
                    print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Measuring responses | Model: {model} | Trained on: {dbt} | Responses to: {dbr} | Image #: {i + 1}/{len(imageFiles)} | Image name: {os.path.basename(outFile)[:-4]}')
                    saveOutputs(model, dbt, paramsFile, imageFile, outFile)
                lastOutFile = outFile

        # run PCA
        for dbc in databases['components']:
            responseDirComp = os.path.join(responseDir, dbc)
            responseFiles = sorted(glob.glob(os.path.join(responseDirComp, '*.pkl')))

            for l, layer in enumerate(layers):
                if layer in theseLayers:

                    outDir = os.path.join('data/DNN', model, dbt, 'PCA', dbc, layer)
                    os.makedirs(outDir, exist_ok=True)

                    pcaFile = os.path.join(outDir, 'pca.pkl')
                    if not os.path.isfile(pcaFile):

                        # collect responses for this layer
                        responses = []
                        for r, responseFile in enumerate(responseFiles):
                            print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Collating responses for PCA | Model: {model} | Layer: {layer} | Trained on: {dbt} | PCA on: {dbc} | Image: {r+1}/{len(responseFiles)}')
                            response = pickle.load(open(responseFile, 'rb'))
                            responses.append(np.array(torch.Tensor.cpu(response[layer].flatten())))
                        #pickle.dump(responses, open(os.path.join(outDir, 'responses.pkl'), 'wb'), protocol=4)
                        # run PCA
                        print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Running PCA | Model: {model} | Layer: {layer} | Trained on: {dbt} | PCA on: {dbc}')

                        pca = PCA()
                        pca.fit(responses)
                        pickle.dump(pca, open(pcaFile, 'wb'), protocol=4) # pickle will fail for filesize over 4GB unless protocol == 4

                        # plot variance explained
                        variance = pca.explained_variance_ratio_
                        cumSumVariance = np.cumsum(variance)
                        plt.ylabel('Proportion of Variance Explained')
                        plt.xlabel('Number of Features')
                        plt.title('PCA Analysis')
                        plt.ylim(0, 1.005)
                        plt.style.context('seaborn-whitegrid')
                        plt.plot(cumSumVariance)
                        plt.savefig(os.path.join(outDir, 'PCAvarExp.png'))
                        plt.show()
                        plt.close()

        # load PCA and get images with highest and lowest weights
        for dbc in databases['components']:
            for l, layer in enumerate(layers):
                if layer in theseLayers:

                    pcaDir = os.path.join('data/DNN', model, dbt, 'PCA', dbc, layer)
                    pcaFile = os.path.join(pcaDir, 'pca.pkl')
                    pca = pickle.load(open(pcaFile, 'rb'))

                    for dbs in databases['samples']:

                        # to leave off at last database analysed, check if last file has been written already
                        lastFile = os.path.join(pcaDir, dbs, f'PC_{nPCs:02}', 'low/tiled.jpg')
                        if not os.path.isfile(lastFile) or overwrite:

                            responseDirSamples = os.path.join(responseDir, dbs)
                            responseFiles = sorted(glob.glob(os.path.join(responseDirSamples, '*.pkl')))
                            imageDir = os.path.join('/home/dave/Datasets', dbs)
                            if not os.path.isdir(imageDir):
                                imageDir = os.path.join('/media/dave/HotProjects/Datasets', dbs)
                            imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.png'), recursive=True))
                            if len(imageFiles) == 0:
                                imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.jpg'), recursive=True))
                            if len(imageFiles) == 0:
                                imageFiles = sorted(glob.glob(os.path.join(imageDir, '**/*.tif*'), recursive=True))
                            if len(imageFiles) == 0:
                                raise Exception('No png, jpg or tif image files found!')
                            nImages = len(imageFiles)

                            # Calculate PC weights for sample images
                            PCweights = np.zeros((nImages, nPCs))
                            imageList = []
                            for r, responseFile in enumerate(responseFiles):
                                print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Calculating PC weights | Model: {model} | Layer: {layer} | Trained on: {dbt} | PCA on: {dbc} | Plot using: {dbs} | Image: {r+1}/{len(responseFiles)}')
                                responseData = pickle.load(open(responseFile, 'rb'))
                                response = np.array([np.array(torch.Tensor.cpu(responseData[layer].flatten()))]) # may need to change between [l] and [l+1] if error is here
                                PCweights[r,:] = pca.transform(response)[:,:nPCs]
                                imageList.append(imageFiles[r])


                            # save copies, make plots for top nPCs
                            print('Selecting images with high magnitude weights and plotting PC weights...')
                            for pc in range(nPCs):
                                rankOrder = PCweights[:, pc].argsort()[::-1]
                                topImageIdxs = rankOrder[:nPCAimages]
                                bottomImageIdxs = rankOrder[-nPCAimages:]
                                for i in range(nPCAimages):

                                    # highest weights
                                    outDir = os.path.join(pcaDir, dbs, f'PC_{pc+1:02}', 'high')
                                    os.makedirs(outDir, exist_ok=True)
                                    topImagePath = imageList[topImageIdxs[i]]
                                    topImageName = os.path.basename(topImagePath[:-4])
                                    outFile = os.path.join(outDir, f'{i+1:02}_{topImageName}.png')
                                    shutil.copy(topImagePath, outFile)
                                    centreCropResize(outFile, outFile, imageSize=[256,256])

                                    # lowest weights
                                    outDir = os.path.join(pcaDir, dbs, f'PC_{pc+1:02}', 'low')
                                    os.makedirs(outDir, exist_ok=True)
                                    bottomImagePath = imageList[bottomImageIdxs[i]]
                                    bottomImageName = os.path.basename(bottomImagePath[:-4])
                                    outFile = os.path.join(outDir, f'{nImages-i:02}_{bottomImageName}.png')
                                    shutil.copy(bottomImagePath, outFile)
                                    centreCropResize(outFile, outFile, imageSize=[256, 256])

                                # make tiled image of images from this PC
                                for highlow in ['high', 'low']:
                                    outDir = os.path.join(pcaDir, dbs, f'PC_{pc+1:02}', highlow)
                                    outImage = os.path.join(outDir, 'tiled.jpg')
                                    images = sorted(glob.glob(os.path.join(outDir, '*.png')))[:nPCAimages]  # list their names
                                    testImage = Image.open(images[0])
                                    imsize = testImage.size
                                    Ncols = int(np.sqrt(len(images)))
                                    Nrows = int(np.sqrt(len(images)))
                                    montage = Image.new(mode='RGB', size=(Nrows * imsize[0], Ncols * imsize[1]), color=0)
                                    for i, im in enumerate(images):
                                        offset_x = (i % Ncols) * imsize[0]
                                        offset_y = int(i / Ncols) * imsize[1]
                                        image = Image.open(im)  # open it
                                        montage.paste(image, (offset_x, offset_y))
                                    montage.save(outImage)

                                    # make another copy in dir for all tiled images
                                    tiledDir = os.path.join(os.path.dirname(pcaDir), 'tiledImages', dbs)
                                    os.makedirs(tiledDir, exist_ok=True)
                                    outImage = os.path.join(tiledDir, f'{layer}_PC{pc+1:02}_{highlow}.jpg')
                                    montage.save(outImage)

                                # make histogram of PCA values across images
                                outDir = os.path.join(pcaDir, dbs, f'PC_{pc + 1:02}')
                                histData = PCweights[:, pc]
                                min = np.min(histData)
                                max = np.max(histData)
                                x_pos = np.linspace(min, max, len(histData))
                                fig, ax = plt.subplots()
                                ax.hist(histData, bins = 32)
                                ax.set_title(f'PC {pc+1}')
                                plt.xlabel('PC weight')
                                plt.ylabel('image count')
                                ax.yaxis.grid(True)
                                plt.tight_layout()
                                plt.savefig(f'{outDir}/PCAvals_hist.png')
                                plt.show()
                                plt.close()

                            # select images from 4 quadrants of first 2 PCS
                            if dbs == 'imagenet1000_subset':
                                nImagesQuad = 25
                                rangePC1 = [np.min(PCweights[:, 0]), np.max(PCweights[:, 0])]
                                rangePC2 = [np.min(PCweights[:, 1]), np.max(PCweights[:, 1])]
                                centreCoords = {'lowlow': [.5 * rangePC1[0], .5 * rangePC2[0]],
                                                'lowhigh': [.5 * rangePC1[0], .5 * rangePC2[1]],
                                                'highlow': [.5 * rangePC1[1], .5 * rangePC2[0]],
                                                'highhigh': [.5 * rangePC1[1], .5 * rangePC2[1]]}
                                for quadrant in centreCoords.keys():
                                    outDir = os.path.join(pcaDir, dbs, 'PC1PC2quadrants', quadrant)
                                    os.makedirs(outDir, exist_ok=True)
                                    distances = np.zeros(nImages)
                                    for i in range(nImages):
                                        distances[i] = np.sqrt(np.sum((centreCoords[quadrant] - PCweights[i,0:2])**2))
                                    rankOrder = distances.argsort()[::-1]
                                    closestImageIdxs = rankOrder[-nImagesQuad:]
                                    for i, im in enumerate(images):
                                        closestImagePath = imageList[closestImageIdxs[i]]
                                        closestImageName = os.path.basename(closestImagePath)
                                        outFile = os.path.join(outDir, f'{i:06}_{closestImageName}')
                                        shutil.copy(closestImagePath, outFile)
                                        centreCropResize(outFile, outFile, imageSize=[256, 256])

                                    # make tiled image of images from this PC
                                    outImage = os.path.join(outDir, f'{quadrant}_tiled.jpg')
                                    images = sorted(glob.glob(os.path.join(outDir, '*.???')))[:nImagesQuad]  # list their names
                                    testImage = Image.open(images[0])
                                    imsize = testImage.size
                                    Ncols = int(np.sqrt(len(images)))
                                    Nrows = int(np.sqrt(len(images)))
                                    montage = Image.new(mode='RGB', size=(Nrows * imsize[0], Ncols * imsize[1]), color=0)
                                    for i, im in enumerate(images):
                                        offset_x = (i % Ncols) * imsize[0]
                                        offset_y = int(i / Ncols) * imsize[1]
                                        image = Image.open(im)  # open it
                                        montage.paste(image, (offset_x, offset_y))
                                    montage.save(outImage)

                            # place BaoTsao3 images in this space
                            if dbs == 'BaoTsao3' and layer == 'fc6':
                                outDir = os.path.join(pcaDir, dbs)
                                PCweights = np.zeros((nImages, nPCs))
                                imageList = []
                                for r, responseFile in enumerate(responseFiles):
                                    print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Calculating PC weights | Model: {model} | Layer: {layer} | Trained on: {dbt} | PCA on: {dbc} | Plot using: {dbs} | Image: {r + 1}/{len(responseFiles)}')
                                    responseData = pickle.load(open(responseFile, 'rb'))
                                    response = np.array([np.array(torch.Tensor.cpu(responseData[layer].flatten()))])
                                    PCweights[r, :] = pca.transform(response)[:, :nPCs]
                                    imageList.append(imageFiles[r])

                                # plot each image
                                PCweights[:,1] = -PCweights[:,1]  # invert so space matches Bao and Tsao
                                fig, ax = plt.subplots(figsize=(3.5, 3.5))
                                conditions = ['spiky-animate', 'stubby-animate', 'spiky-inanimate', 'stubby-inanimate']
                                colours = ['chartreuse', 'deepskyblue', 'orange','fuchsia']
                                for c in range(4):
                                    imageIdx = np.arange((c*20),(c*20+20))
                                    ax.scatter(PCweights[imageIdx, 0], PCweights[imageIdx, 1], s=30, color=colours[c])

                                # plot mean of points
                                for c in range(4):
                                    imageIdx = np.arange((c * 20), (c * 20 + 20))
                                    mean1 = np.mean(PCweights[imageIdx, 0])
                                    mean2 = np.mean(PCweights[imageIdx, 1])
                                    ax.scatter(mean1, mean2, s=100, edgecolors='black', color=colours[c], label=conditions[c])

                                ax.set_title(f'object space stimuli\nin PC1-PC2 space')
                                plt.xlabel('PC1')
                                plt.ylabel('PC2')
                                plt.axvline(0, c='black')
                                plt.axhline(0, c='black')
                                plt.xlim(-70, 70)
                                plt.ylim(-70, 70)
                                plt.tight_layout()
                                #plt.legend(frameon=False)
                                plt.savefig(f'{outDir}/AARstimuli.png')
                                plt.show(block=False)
                                plt.close()

                            # place BFHO images in this space
                            if dbs == 'BFHO' and layer == 'fc6':
                                outDir = os.path.join(pcaDir, dbs)
                                PCweights = np.zeros((nImages, nPCs))
                                imageList = []
                                for r, responseFile in enumerate(responseFiles):
                                    print(
                                        f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Calculating PC weights | Model: {model} | Layer: {layer} | Trained on: {dbt} | PCA on: {dbc} | Plot using: {dbs} | Image: {r + 1}/{len(responseFiles)}')
                                    responseData = pickle.load(open(responseFile, 'rb'))
                                    response = np.array([np.array(torch.Tensor.cpu(responseData[layer].flatten()))])
                                    PCweights[r, :] = pca.transform(response)[:, :nPCs]
                                    imageList.append(imageFiles[r])

                                # plot each image
                                PCweights[:,1] = -PCweights[:,1]  # invert so space matches Bao and Tsao
                                fig, ax = plt.subplots(figsize=(3.5, 3.5))
                                conditions = ['body', 'face', 'house', 'object']
                                colours = ['chartreuse', 'deepskyblue', 'red', 'gold']
                                for c, cond in enumerate(conditions):
                                    theseImagePaths = []
                                    imageIdx = []
                                    for i, image in enumerate(imageList):
                                        if cond in image:
                                            theseImagePaths.append(image)
                                            imageIdx.append(i)
                                    ax.scatter(PCweights[imageIdx, 0], PCweights[imageIdx, 1], s=30, color=colours[c])

                                # plot mean of points
                                for c, cond in enumerate(conditions):
                                    theseImagePaths = []
                                    imageIdx = []
                                    for i, image in enumerate(imageList):
                                        if cond in image:
                                            theseImagePaths.append(image)
                                            imageIdx.append(i)
                                    mean1 = np.mean(PCweights[imageIdx, 0])
                                    mean2 = np.mean(PCweights[imageIdx, 1])
                                    ax.scatter(mean1, mean2, s=100, edgecolors='black', color=colours[c], label=cond)

                                ax.set_title(f'object categories\nin PC1-PC2 space')
                                plt.xlabel('PC1')
                                plt.ylabel('PC2')
                                plt.axvline(0, c='black')
                                plt.axhline(0, c='black')
                                plt.xlim(-70, 70)
                                plt.ylim(-70, 70)
                                plt.tight_layout()
                                #plt.legend(frameon=False)
                                plt.savefig(f'{outDir}/BFHOstimuli.png')
                                plt.show(block=False)
                                plt.close()




        # delete responses to save a lot of disk space
        #shutil.rmtree(responseDir)



