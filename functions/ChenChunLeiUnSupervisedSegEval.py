# -*- coding: cp936 -*-
import os, re, sys, math, arcpy, numpy as np
from array import array
from trim_fea_class_names import trim_feature_class
from sort_feature_names import FeatureQuickSort
from get_fea_names import export_fea_names

curDir = os.getcwd()
arcpy.env.workspace = curDir + '\\results'
arcpy.env.overwriteOutput=True

TOTAL_AREA = 2.5e5
FIRSTIDENTIFY_NUM = 0

feature_names = export_fea_names(curDir)
feature_names,HEADERIDENTIFY_NUM = trim_feature_class(feature_names,[1,1,127],shp=[1,3,5,7,9],cpt =[1,3,5,7,9],filename = 'ChenChunLeiSegEvalDoneFeas.txt')
FeatureQuickSort(feature_names)

outFeatureClass = curDir + '\\auxiliary_files.gdb\\featurePolygon'
pattern = re.compile(r'\d+')
feaNameFields = [u'FID', u'Shape', u'Mean_1', u'Mean_2', u'Mean_3', u'Mean_4', u'Mean_5', u'Mean_6', u'Mean_7', \
    u'Mean_8', u'Standard_1', u'Standard_2', u'Standard_3', u'Standard_4',\
     u'Standard_5', u'Standard_6', u'Standard_7', u'Standard_8',\
     u'SpectVar', u'ShpPeri', u'ShpArea', u'PixelNum', u'FIDOrigin']
calFieldlist = [ element for element in feaNameFields if pattern.findall(element) or u'FID' == element]        
calFieldlist.append(u'SHAPE@LENGTH'); calFieldlist.append(u'SHAPE@AREA')
calFieldlistInner = [ element for element in feaNameFields if u'Mean' in element or u'FID' == element]   
print feature_names
for featureName in feature_names:
    arcpy.PolygonToLine_management(featureName, outFeatureClass)
    #  segmentation index 
    ASEI = 0
    cursor = arcpy.da.SearchCursor(featureName,calFieldlist,sql_clause = (None, 'ORDER BY FID'))
    for row in cursor:
        # line length
        sql_expression = 'LEFT_FID = {0} or RIGHT_FID = {1}'.format(row[0], row[0])
        arcpy.MakeFeatureLayer_management (outFeatureClass, "FeaturePolyline")
        arcpy.SelectLayerByAttribute_management ("FeaturePolyline", "NEW_SELECTION", sql_expression)
        LineFIDSet = set()
        LineLengthDict = {}
        
        loc_cursorline = arcpy.da.SearchCursor("FeaturePolyline",['LEFT_FID', 'RIGHT_FID', 'SHAPE@LENGTH'])
        for loc_row in loc_cursorline:
            LineFIDSet = LineFIDSet.union(set(loc_row[0:2])) - {-1} - {row[0]}
            LineFID = list(set(loc_row[0:2]) - {row[0]})[0]
            LineLengthDict[LineFID] = loc_row[2]/row[calFieldlist.index(u'SHAPE@LENGTH')]
        LineFIDSet = list(LineFIDSet)
        
        # area
        arcpy.MakeFeatureLayer_management (featureName, "FeaturePolygon")
        for i in range(len(LineFIDSet)):
            sql_expression = 'FID = {0}'.format(LineFIDSet[i])
            arcpy.SelectLayerByAttribute_management ("FeaturePolygon", "ADD_TO_SELECTION", sql_expression)
        loc_cursorpoly = arcpy.da.SearchCursor("FeaturePolygon", calFieldlistInner)
        
        heteroIndex = np.zeros([1,3])
        for loc_row in loc_cursorpoly:
            meanValue = array("d")
            meanValue = loc_row[calFieldlistInner.index(u'Mean_1'):(calFieldlistInner.index(u'Mean_3')+1)]
            meanValue = np.array(meanValue, dtype = np.float)
            meanValue = meanValue * LineLengthDict[loc_row[0]]
            heteroIndex = np.concatenate([heteroIndex,[meanValue]])
        heteroIndex = np.delete(heteroIndex, 0, 0)  
        heteroIndex = np.sum(heteroIndex)
        homoIndex = array("d")
        [homoIndex.append(row[i]) for i in range(calFieldlist.index(u'Standard_1'),  calFieldlist.index(u'Standard_3')+1)]
        homoIndex = np.frombuffer(homoIndex, dtype = np.float)
        # weight
        weight = np.array([1,1,1])
        # merge the eight layers
        SEI = np.sum(heteroIndex / homoIndex * weight)
        ASEI = ASEI + row[calFieldlist.index(u'SHAPE@AREA')] * SEI
    
    featurename = featureName.split('.')[0]
    print featurename
    OutputFileName = curDir + '\\evaluation_results\\'+'ChenChunLeiUnsuperviseSegEval.csv'
    fileID = open(OutputFileName,'a+')
    if FIRSTIDENTIFY_NUM == 0 and HEADERIDENTIFY_NUM ==0:
        fileID.write('featureName,ChenChunLeiUnsuperviseSegEval')
        fileID.write('\n')
        fileID.write('{0:^26},{1:9.6g}'.format(featurename,ASEI/TOTAL_AREA))
        fileID.write('\n')
    else:
        fileID.write('{0:^26},{1:9.6g}'.format(featurename,ASEI/TOTAL_AREA))
        fileID.write('\n')
    fileID.close()
    FIRSTIDENTIFY_NUM = FIRSTIDENTIFY_NUM + 1
    
    
    file_feature_names = open(curDir + '\\evaluation_results\\UnSupervisedSegEvalFeas.txt', 'a+')
    file_feature_names.write(featurename)
    file_feature_names.write('\n')
    file_feature_names.close()