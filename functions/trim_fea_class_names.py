import os 
import sys
import numpy as np
import re

def trim_feature_class(feature_names, scale = [1,1,127], shp = [1, 3, 5, 7, 9], cpt = [1, 3, 5, 7, 9], filename = 'DoneSupervisedSegEvalFeas.txt', spect = [1,5,10,15,20,25]):
    """ 
        full path of the done file names
    """
    p = re.compile('\d+')
    scale = np.arange(scale[0],scale[2],scale[1]).tolist()
    featureNameList = []
    for featureName in feature_names:
        temp = p.findall(featureName)
        if int(temp[0]) in shp and int(temp[1]) in cpt and int(temp[2]) in scale:
            featureNameList.append(featureName)
    filenameID = open(filename, 'a+')
    fileNameList = []
    for fileStr in filenameID.readlines():
        fileNameList.append(fileStr.split('\n')[0] + '.shp')
    featureNameSet = set(featureNameList) - set(fileNameList)
    featureNameList = list(featureNameSet)
    return featureNameList, len(fileNameList)
    
