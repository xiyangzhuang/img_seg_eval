import cPickle as pickle
import os, sys, re, time, glob, win32api
from functions.trim_fea_class_names import trim_feature_class
from functions.sort_feature_names import FeatureQuickSort
from functions.segassess.SegmentationAssess0424AM import segment_evaluation, thesis_segment_evaluation, rand_sample
from functions.get_fea_names import export_fea_names
# global variable to use
namepattern = re.compile('\dand\dScale\d+', re.I)

# mark the file generate, code block generate file stamp
locTime = time.localtime()
segFileList = glob.glob('functions\\segassess\\SegmentationAssess*.py')
segFileList.sort()
pattern = re.compile('\d{4,8}[AP]\w')
strStamp = pattern.findall(segFileList.pop())[0]

curDir = os.path.dirname(os.path.dirname(__file__))
resultDir = curDir + '\\' + 'evaluation_results\\R_' + strStamp
if not os.path.isdir(resultDir):
    os.mkdir(resultDir)
else:
    if len(os.listdir(resultDir)) <> 0:
        strStamp = strStamp + '0' + str(locTime.tm_mon) + str(locTime.tm_mday) + str(locTime.tm_hour)

def exitWithState(sig, func=None):
    OutAFIFile = file(resultDir + '\\AFI{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictAFI, OutAFIFile)
    OutAFIFile.close()
    OutOLFile = file(resultDir + '\\OL{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictOL, OutOLFile)
    OutOLFile.close()
    OutIFile = file(resultDir + '\\I{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictI, OutIFile)    
    OutIFile.close()
    raise KeyboardInterrupt
win32api.SetConsoleCtrlHandler(exitWithState, 1)

#get the feature names need to be calculated
feature_names = export_fea_names(curDir)
feature_names,HEADERIDENTIFY_NUM = trim_feature_class(feature_names,[1,1,127],shp=[1,3,5,7,9],cpt =[1,3,5,7,9],
                                    filename = resultDir + '\\Calculated{0}{1}{2}.txt'.format(max(arrObjectID),className[0],strStamp))
#feature_names,HEADERIDENTIFY_NUM = trim_feature_class(feature_names,[10,10,127],shp=[1],cpt =[5],filename = 'CalculatedFeas.txt')
FeatureQuickSort(feature_names)

FIRSTIDENTIFY_NUM = 0
dictAFI = {}; dictOL = {}; dictI = {};
try:
    for featureName in feature_names:
        intersect_area = 0
        # AFI, OL, I reversion export
        arrObjectID = range(1,34); className = ['all', 'building', 'tree', 'water']
        OutputFileName = resultDir + '\\R_{0}{1}{2}.csv'.format(max(arrObjectID),className[0],strStamp)
        dictAttr = rand_sample(objectIDList = arrObjectID, class_name = className, sample_rate = 1)
        OverallUSE, OverallOSE, OverallBDI, OverallPDI, OverallAFI, OverallOL, OverallI, AFI, OL, I = thesis_segment_evaluation(featureName, intersect_area, *dictAttr)
        featurename = featureName.split('.')[0]
        if os.path.isfile(resultDir + '\\AFI{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp)) and FIRSTIDENTIFY_NUM==0:
            dictAFIFile = open(resultDir + '\\AFI{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'rb')
            dictOLFile = open(resultDir + '\\OL{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'rb')
            dictIFile = open(resultDir + '\\I{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'rb')
            dictAFI = pickle.load(dictAFIFile)
            dictOL = pickle.load(dictOLFile)
            dictI = pickle.load(dictIFile)
        else:
            keyname = namepattern.findall(featurename)[0]
            dictAFI[keyname] = AFI
            dictOL[keyname] = OL
            dictI[keyname] = I
        fileID = open(OutputFileName,'a+')
        if HEADERIDENTIFY_NUM ==0 and FIRSTIDENTIFY_NUM==0:
            fileID.write('featureName,OverallUSE, OverallOSE, OverallBDI, OverallPDI, OverallAFI, OverallOL, OverallI')
            fileID.write('\n')
            fileID.write('{0:^26},{1:9.6g},{2:9.6g},{3:9.6g},{4:9.6g},{5:9.6g},{6:9.6g},{7:9.6g}'.format(featurename,OverallUSE, OverallOSE, OverallBDI, OverallPDI, OverallAFI, OverallOL, OverallI))
            fileID.write('\n')
        else:
            fileID.write('{0:^26},{1:9.6g},{2:9.6g},{3:9.6g},{4:9.6g},{5:9.6g},{6:9.6g},{7:9.6g}'.format(featurename,OverallUSE, OverallOSE, OverallBDI, OverallPDI, OverallAFI, OverallOL, OverallI))
            fileID.write('\n')
        fileID.close()
        FIRSTIDENTIFY_NUM = FIRSTIDENTIFY_NUM + 1
        ### if the file name is not in the range, it means it calculate all the reference polygon
        file_feature_names = open(resultDir + '\\Calculated{0}{1}{2}.txt'.format(max(arrObjectID),className[0],strStamp), 'a+')
        file_feature_names.write(featurename)
        file_feature_names.write('\n')
        file_feature_names.close()
    OutAFIFile = file(resultDir + '\\AFI{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictAFI, OutAFIFile)
    OutAFIFile.close()
    OutOLFile = file(resultDir + '\\OL{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictOL, OutOLFile)
    OutOLFile.close()
    OutIFile = file(resultDir + '\\I{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictI, OutIFile)    
    OutIFile.close()
  
except KeyboardInterrupt:
    OutAFIFile = file(resultDir + '\\AFI{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictAFI, OutAFIFile)
    OutAFIFile.close()
    OutOLFile = file(resultDir + '\\OL{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictOL, OutOLFile)
    OutOLFile.close()
    OutIFile = file(resultDir + '\\I{0}{1}{2}.pkl'.format(max(arrObjectID),className[0],strStamp), 'wb')
    pickle.dump(dictI, OutIFile)    
    OutIFile.close()
#exec(open(curDir + "\\" + "MainFunctionDiscrepancyIndex.py").read())
    
