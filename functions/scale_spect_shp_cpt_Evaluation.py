# -*- coding: cp936 -*-
import arcpy, os, re, sys, math, random
import numpy as np
from array import array
from trim_fea_class_names import trim_feature_class

curDir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
arcpy.env.workspace = curDir + '\\results'
arcpy.env.overwriteOutput=True
workspace = curDir + '\\results'
TOTAL_AREA = 2.5e5
TOTAL_PIXEL_NUM = 1e6


rPoly = curDir + '\\reference_features.gdb\\ReferencyPoly'
rPointSelected = curDir + '\\auxiliary_files.gdb\\ReferencyPointRateSelected'
rCentroid = curDir + '\\reference_features.gdb\\ReferencyPolyPoint'

spatial_ref = arcpy.Describe(rCentroid).spatialReference
                                
def slim_feature(featureName, listReferPointAttr):
    arcpy.MakeFeatureLayer_management(featureName, "feature_layer")
    arcpy.SelectLayerByLocation_management("feature_layer","INTERSECT","rPoly","","NEW_SELECTION")
    arcpy.Identity_analysis("rPoly","feature_layer",identify_feature)
    
def comborand_sample(sample_rate = 1, class_name = ['tree', 'building', 'water', 'car', 'other']):
    numDIR = 50
    sampleNum = int(sample_rate * numDIR)
    sampleID = random.sample(range(1,numDIR+1), sampleNum)
    dictReferPolyAttr = {row[0]:[row[1],row[2]] for row in arcpy.da.SearchCursor(rPoly,["OBJECTID","SHAPE_Area","SHAPE@","CLASSNAME"]) if row[0] in sampleID and row[3] in class_name};
    dictReferPointAttr = {row[0]:[row[1],row[2]] for row in arcpy.da.SearchCursor(rCentroid,["OBJECTID","ORIG_FID","SHAPE@","CLASSNAME"]) if row[3] in class_name and row[0] in sampleID};
    listReferPointAttr = [row[2] for row in arcpy.da.SearchCursor(rCentroid,["OBJECTID","ORIG_FID","SHAPE@","CLASSNAME"]) if row[3] in class_name and row[0] in sampleID];
    arcpy.CopyFeatures_management(listReferPointAttr, rPointSelected)
    ReferPolyWeights = np.ones(len(dictReferPolyAttr))
    return dictReferPolyAttr, dictReferPointAttr, ReferPolyWeights
    
def normalize(iterable_object, method_type = 0):
    normResults = array("d")
    if method_type == 0:
        if np.std(iterable_object) <> 0:
            [normResults.append((i-np.average(iterable_object)) / np.std(iterable_object)) for i in iterable_object if not np.isnan(i)]
        else:
            [normResults.append(i-np.average(iterable_object)) for i in iterable_object if not np.isnan(i)]
    elif method_type == 1:
        minValue = np.nanmin(iterable_object)
        maxValue = np.nanmax(iterable_object)
        if maxValue == minValue:
            [normResults.append(i - minValue) for i in iterable_object if not np.isnan(i)]
        else:
           [normResults.append((i - minValue)/(maxValue - minValue)) for i in iterable_object if not np.isnan(i)] 
        
    return normResults 
    
def multiweight(iterable_object, weights):
    multiweightResults = array("d")
    total = sum(weights)
    if len(iterable_object) == 0:
        multiweightResults.append(0)
        return multiweightResults
    else:
        [multiweightResults.append(i*j/total) for i in iterable_object for j in weights]
        return multiweightResults   
        
def combosegment_evaluation(feature_name, intersect_area, dictReferPolyAttr, dictReferPointAttr, ReferPolyWeights):
    """
        dictRefer include two dictionary element: dictReferPolyAttr,dictReferPointAttr
    """
    # data process
    identify_feature = curDir + '\\auxiliary_files.gdb\\identifiedFeature'
    identifySelect_feature = curDir + '\\auxiliary_files.gdb\\identifiedExtractFeature'
    arcpy.MakeFeatureLayer_management(rPoly, "rPoly")
    arcpy.SelectLayerByLocation_management("rPoly","INTERSECT",rPointSelected,"","NEW_SELECTION")
    arcpy.MakeFeatureLayer_management(feature_name, "feature_layer")
    arcpy.SelectLayerByLocation_management("feature_layer","INTERSECT","rPoly","","NEW_SELECTION")
    dictFeaAttr = {row[0]:[row[1],row[2]] for row in arcpy.da.SearchCursor("feature_layer",["FID","SHAPE@AREA","SHAPE@"])}
    arcpy.Identity_analysis("rPoly","feature_layer",identify_feature)
    
    
    pattern = re.compile('^FID_')
    SortedField = [field.name for field in arcpy.ListFields(identify_feature) if pattern.match(field.name)]
    CalFields = [field.name for field in arcpy.ListFields(identify_feature)]
    indAIR = CalFields.index(SortedField[1])
    indDIR = CalFields.index(SortedField[0])
    indArea = CalFields.index('SHAPE_Area')
    indLength = CalFields.index('SHAPE_Length')
    # select the feature coincide with the intersect_area we set
    if intersect_area <> 0:
        listIdentityAttr = [row[2] for row in arcpy.da.SearchCursor(identify_feature,[ SortedField[0] ,"SHAPE@AREA","SHAPE@"]) if row[1]/dictReferPolyAttr[row[0]][0] >= intersect_area]
        listfeaindDIR = [row[0] for row in arcpy.da.SearchCursor(identify_feature,[SortedField[0] ,"SHAPE@AREA","SHAPE@"]) if row[1]/dictReferPolyAttr[row[0]][0] >= intersect_area]
        arcpy.CopyFeatures_management(listIdentityAttr, identifySelect_feature)
        arcpy.MakeFeatureLayer_management(identify_feature, 'identify_feature')
        arcpy.SelectLayerByLocation_management('identify_feature',"ARE_IDENTICAL_TO",identifySelect_feature,"","NEW_SELECTION")
        arcpy.SelectLayerByLocation_management('feature_layer',"INTERSECT",identifySelect_feature,"","NEW_SELECTION")
        cursor = arcpy.da.SearchCursor('identify_feature', CalFields, sql_clause=(None, 'ORDER BY {0} DESC'.format(SortedField[0])))
    else:
        cursor = arcpy.da.SearchCursor(identify_feature, CalFields, sql_clause=(None, 'ORDER BY {0} DESC'.format(SortedField[0])))
        
    
    row = cursor.next(); initFID = row[indDIR]; 
    cursor.reset()
    intersectFeaArea = {}
    underMergingIJ = array("d"); overMergingIJ = array("d"); underMergingI = array("d"); overMergingI = array("d");
    underMergingIW = array("d"); overMergingIW = array("d");
    
    simsizeIJ = array("d"); simsizeI = array("d"); simsizeIW = array("d"); 
    
    rasubIJ = array("d"); rasuperIJ = array("d"); rasubI = array("d"); rasuperI = array("d"); rasubIW = array("d"); rasuperIW = array("d"); rasubIJW = array("d"); rasuperIJW = array("d")
    
    qrIJ = array("d"); qrI = array("d"); qrIW = array("d");
    
    oversegmentationIJ = array("d"); undersegmentationIJ = array("d"); oversegmentationIJW = array("d"); undersegmentationIJW = array("d"); 
    oversegmentationI = array("d"); undersegmentationI = array("d"); oversegmentationIW = array("d"); undersegmentationIW = array("d")
    
    qlocIJ = array("d"); rpsubIJ = array("d"); rpsuperIJ = array("d")
    qlocI = array("d"); rpsubI = array("d"); rpsuperI = array("d"); qlocIW = array("d"); rpsubIW = array("d"); rpsuperIW = array("d")
    
    afiI = array("d"); countoverI = array("d"); countunderI = array("d")
    # combo method
    mergesumI = array("d")
    overunderI = array("d");
    dIJ = array("d"); dIJW = array("d");  mergesumIJW = array("d"); dIW = array("d"); dI = array("d")
    
    MIJ = array("d"); MIJW = array("d"); MIW = array("d"); MI = array("d")
    
    zh1IJ = array("d"); zh1IJW = array("d"); zh1IW = array("d"); zh1I = array("d") 
    zh2IJ = array("d"); zh2IJW = array("d"); zh2IW = array("d"); zh2I = array("d")
    for row in cursor:
        if row[indDIR] == initFID:
            underMergingIJ.append( (dictReferPolyAttr[row[indDIR]][0] - row[indArea])/ dictReferPolyAttr[row[indDIR]][0] ); 
            overMergingIJ.append( (dictFeaAttr[row[indAIR]][0] - row[indArea])/ dictReferPolyAttr[row[indDIR]][0] );
            intersectFeaArea[dictFeaAttr[row[indAIR]][0]] = row[indAIR]
            simsizeIJ.append(min(dictReferPolyAttr[row[indDIR]][0], dictFeaAttr[row[indAIR]][0]) / max(dictReferPolyAttr[row[indDIR]][0], dictFeaAttr[row[indAIR]][0]))
            rasubIJ.append(row[indArea]/dictReferPolyAttr[row[indDIR]][0]); rasuperIJ.append(row[indArea]/dictFeaAttr[row[indAIR]][0])
            qrIJ.append(1- row[indArea]/(dictReferPolyAttr[row[indDIR]][0]+dictFeaAttr[row[indAIR]][0]-row[indArea]))
            oversegmentationIJ.append(1- row[indArea]/dictReferPolyAttr[row[indDIR]][0]);undersegmentationIJ.append(1- row[indArea]/dictFeaAttr[row[indAIR]][0])
            distance = dictReferPointAttr[row[indDIR]][1].distanceTo(arcpy.PointGeometry(
                        dictFeaAttr[row[indAIR]][1].trueCentroid, spatial_ref) )
            qlocIJ.append(distance)
            rpsubIJ = qlocIJ;
            
        elif row[indDIR] <> initFID:
            if 'unweighted':
                rpsuperIJ = np.divide(rpsubIJ, np.max(rpsubIJ))
                
                underMergingI.append(np.average(underMergingIJ));overMergingI.append(np.average(overMergingIJ))
                simsizeI.append(np.average(simsizeIJ)); 
                rasubI.append(np.average(rasubIJ)); rasuperI.append(np.average(rasuperIJ))
                qrI.append(np.average(qrIJ))
                oversegmentationI.append(np.average(oversegmentationIJ)); undersegmentationI.append(np.average(undersegmentationIJ))
                qlocI.append(np.average(qlocIJ)); rpsubI.append(np.average(rpsubIJ)); rpsuperI.append(np.average(rpsuperIJ))
                # not weighted combo index 
                dIJ = np.sqrt(np.divide(np.add(np.power(oversegmentationIJ, 2) ,np.power(undersegmentationIJ, 2)), 2))
                MIJ = np.sqrt(np.divide( 
                    np.add(
                        np.add( np.power(np.add(rasubIJ, -1), 2), np.power(np.add(rasuperIJ, -1), 2)), 
                        np.add( np.power(normalize(rpsubIJ), 2), np.power(normalize(rpsuperIJ), 2))
                    ), 4)
                )
                zh1IJ = np.sqrt( np.divide(
                    np.add(
                        np.add( np.power( np.add(simsizeIJ, -1), 2), np.power(normalize(simsizeIJ), 2)),
                        np.add( np.power(normalize(qlocIJ),2), np.power(np.sqrt(qlocIJ), 2))
                    ), 4)
                )
                zh2IJ = np.sqrt( np.divide(np.add(
                    np.power(np.add(simsizeIJ, -1), 2), np.power(normalize(qlocIJ), 2)
                    ), 2)
                )
                dI.append(np.average(dIJ)); MI.append(np.average(MIJ)); zh1I.append(np.average(zh1IJ)); zh2I.append(np.average(zh2IJ));
                
            if 'weighted':
                # combo temp vars
                oversegmentationIJW = multiweight(oversegmentationIJ, weights = qlocIJ); undersegmentationIJW = multiweight(undersegmentationIJ, weights = qlocIJ)
                rasubIJW = multiweight(rasubIJ, weights = qlocIJ); rasuperIJW = multiweight(rasuperIJ, weights = qlocIJ)
                rpsubIJW = multiweight(rpsubIJ, weights = qlocIJ); rpsuperIJW = multiweight(rpsuperIJ, weights = qlocIJ);
                simsizeIJW = multiweight(simsizeIJ, weights = qlocIJ);
                qlocIJW = multiweight(qlocIJ, weights = qlocIJ);
                
                underMergingIW.append(np.average(underMergingIJ, weights = qlocIJ));overMergingIW.append(np.average(overMergingIJ,weights = qlocIJ))
                simsizeIW.append(np.average(simsizeIJ, weights = qlocIJ)); 
                rasubIW.append(np.average(rasubIJ, weights = qlocIJ)); rasuperIW.append(np.average(rasuperIJ, weights = qlocIJ))
                qrIW.append(np.average(qrIJ, weights = qlocIJ))
                oversegmentationIW.append(np.average(oversegmentationIJ, weights = qlocIJ)); undersegmentationIW.append(np.average(undersegmentationIJ, weights = qlocIJ))
                # weighted combo method
                dIJW = np.sqrt(np.divide(np.add(np.power(oversegmentationIJW, 2) ,np.power(undersegmentationIJW, 2)), 2))
                MIJW = np.sqrt(np.divide( 
                    np.add(
                        np.add( np.power(np.add(rasubIJW, -1), 2), np.power(np.add(rasuperIJW, -1), 2)), 
                        np.add( np.power(normalize(rpsubIJW), 2), np.power(normalize(rpsuperIJW), 2))
                    ), 4)
                )
                zh1IJW = np.sqrt( np.divide(
                    np.add(
                        np.add( np.power( np.add(simsizeIJW, -1), 2), np.power(normalize(simsizeIJW), 2)),
                        np.add( np.power(normalize(qlocIJW),2), np.power(np.sqrt(qlocIJW), 2))
                    ), 4)
                )
                zh2IJW = np.sqrt( np.divide(np.add(
                    np.power(np.add(simsizeIJW, -1), 2), np.power(normalize(qlocIJW), 2)
                    ), 2)
                )
                dIW.append(np.average(dIJW)); MIW.append(np.average(MIJW)); zh1IW.append(np.average(zh1IJW)); zh2IW.append(np.average(zh2IJW));
            mergesumI.append(np.average(underMergingIJ, weights = qlocIJ) + np.average(overMergingIJ, weights = qlocIJ))    
            temp_afiI = (dictReferPolyAttr[row[indDIR]][0] - max(intersectFeaArea.keys()))/dictReferPolyAttr[row[indDIR]][0]
            afiI.append(temp_afiI);
            if not dictReferPolyAttr[row[indDIR]][1].contains(dictFeaAttr[intersectFeaArea[max(intersectFeaArea.keys())]][1]) and temp_afiI > 0:
                countunderI.append(1)
            elif dictReferPolyAttr[row[indDIR]][1].contains(dictFeaAttr[intersectFeaArea[max(intersectFeaArea.keys())]][1]) and temp_afiI < 0:
                countoverI.append(1)
            initFID = row[indDIR]
            # initialize the variables
            underMergingIJ = array("d"); overMergingIJ = array("d");simsizeIJ = array("d"); rasubIJ = array("d"); rasuperIJ = array("d");
            qrIJ = array("d");oversegmentationIJ = array("d");undersegmentationIJ = array("d"); qlocIJ = array("d"); rpsubIJ = array("d"); qlocIJ = array("d");
            #### the same code of the first logic ####
            underMergingIJ.append( (dictReferPolyAttr[row[indDIR]][0] - row[indArea])/ dictReferPolyAttr[row[indDIR]][0] ); 
            overMergingIJ.append( (dictFeaAttr[row[indAIR]][0] - row[indArea])/ dictReferPolyAttr[row[indDIR]][0] );
            intersectFeaArea[dictFeaAttr[row[indAIR]][0]] = row[indAIR]
            simsizeIJ.append(min(dictReferPolyAttr[row[indDIR]][0], dictFeaAttr[row[indAIR]][0]) / max(dictReferPolyAttr[row[indDIR]][0], dictFeaAttr[row[indAIR]][0]))
            rasubIJ.append(row[indArea]/dictReferPolyAttr[row[indDIR]][0]); rasuperIJ.append(row[indArea]/dictFeaAttr[row[indAIR]][0])
            qrIJ.append(1- row[indArea]/(dictReferPolyAttr[row[indDIR]][0]+dictFeaAttr[row[indAIR]][0]-row[indArea]))
            oversegmentationIJ.append(1- row[indArea]/dictReferPolyAttr[row[indDIR]][0]);undersegmentationIJ.append(1- row[indArea]/dictFeaAttr[row[indAIR]][0])
            distance = dictReferPointAttr[row[indDIR]][1].distanceTo(arcpy.PointGeometry(
                        dictFeaAttr[row[indAIR]][1].trueCentroid, spatial_ref) )
            qlocIJ.append(distance)
            rpsubIJ = qlocIJ;
            #### the same code of the first logic ####
    #### the same code of the second logic ####
    if 'unweighted':
        rpsuperIJ = np.divide(rpsubIJ, np.max(rpsubIJ))
        underMergingI.append(np.average(underMergingIJ));overMergingI.append(np.average(overMergingIJ))
        simsizeI.append(np.average(simsizeIJ)); 
        rasubI.append(np.average(rasubIJ)); rasuperI.append(np.average(rasuperIJ))
        qrI.append(np.average(qrIJ))
        oversegmentationI.append(np.average(oversegmentationIJ)); undersegmentationI.append(np.average(undersegmentationIJ))
        qlocI.append(np.average(qlocIJ)); rpsubI.append(np.average(rpsubIJ)); rpsuperI.append(np.average(rpsuperIJ))
        # not weighted combo index 
        dIJ = np.sqrt(np.divide(np.add(np.power(oversegmentationIJ, 2) ,np.power(undersegmentationIJ, 2)), 2))
        MIJ = np.sqrt(np.divide( 
            np.add(
                np.add( np.power(np.add(rasubIJ, -1), 2), np.power(np.add(rasuperIJ, -1), 2)), 
                np.add( np.power(normalize(rpsubIJ), 2), np.power(normalize(rpsuperIJ), 2))
            ), 4)
        )
        zh1IJ = np.sqrt( np.divide(
            np.add(
                np.add( np.power( np.add(simsizeIJ, -1), 2), np.power(normalize(simsizeIJ), 2)),
                np.add( np.power(normalize(qlocIJ),2), np.power(np.sqrt(qlocIJ), 2))
            ), 4)
        )
        zh2IJ = np.sqrt( np.divide(np.add(
            np.power(np.add(simsizeIJ, -1), 2), np.power(normalize(qlocIJ), 2)
            ), 2)
        )
        dI.append(np.average(dIJ)); MI.append(np.average(MIJ)); zh1I.append(np.average(zh1IJ)); zh2I.append(np.average(zh2IJ));   
    if 'weighted':
        # combo temp vars
        oversegmentationIJW = multiweight(oversegmentationIJ, weights = qlocIJ); undersegmentationIJW = multiweight(undersegmentationIJ, weights = qlocIJ)
        rasubIJW = multiweight(rasubIJ, weights = qlocIJ); rasuperIJW = multiweight(rasuperIJ, weights = qlocIJ)
        rpsubIJW = multiweight(rpsubIJ, weights = qlocIJ); rpsuperIJW = multiweight(rpsuperIJ, weights = qlocIJ);
        simsizeIJW = multiweight(simsizeIJ, weights = qlocIJ);
        qlocIJW = multiweight(qlocIJ, weights = qlocIJ);
        
        underMergingIW.append(np.average(underMergingIJ, weights = qlocIJ));overMergingIW.append(np.average(overMergingIJ,weights = qlocIJ))
        simsizeIW.append(np.average(simsizeIJ, weights = qlocIJ)); 
        rasubIW.append(np.average(rasubIJ, weights = qlocIJ)); rasuperIW.append(np.average(rasuperIJ, weights = qlocIJ))
        qrIW.append(np.average(qrIJ, weights = qlocIJ))
        oversegmentationIW.append(np.average(oversegmentationIJ, weights = qlocIJ)); undersegmentationIW.append(np.average(undersegmentationIJ, weights = qlocIJ))
        # weighted combo method
        dIJW = np.sqrt(np.divide(np.add(np.power(oversegmentationIJW, 2) ,np.power(undersegmentationIJW, 2)), 2))
        MIJW = np.sqrt(np.divide( 
            np.add(
                np.add( np.power(np.add(rasubIJW, -1), 2), np.power(np.add(rasuperIJW, -1), 2)), 
                np.add( np.power(normalize(rpsubIJW), 2), np.power(normalize(rpsuperIJW), 2))
            ), 4)
        )
        zh1IJW = np.sqrt( np.divide(
            np.add(
                np.add( np.power( np.add(simsizeIJW, -1), 2), np.power(normalize(simsizeIJW), 2)),
                np.add( np.power(normalize(qlocIJW),2), np.power(np.sqrt(qlocIJW), 2))
            ), 4)
        )
        zh2IJW = np.sqrt( np.divide(np.add(
            np.power(np.add(simsizeIJW, -1), 2), np.power(normalize(qlocIJW), 2)
            ), 2)
        )
        dIW.append(np.average(dIJW)); MIW.append(np.average(MIJW)); zh1IW.append(np.average(zh1IJW)); zh2IW.append(np.average(zh2IJW));
    mergesumI.append(np.average(underMergingIJ, weights = qlocIJ) + np.average(overMergingIJ, weights = qlocIJ))      
    temp_afiI = (dictReferPolyAttr[row[indDIR]][0] - max(intersectFeaArea.keys()))/dictReferPolyAttr[row[indDIR]][0]
    afiI.append(temp_afiI);
    if not dictReferPolyAttr[row[indDIR]][1].contains(dictFeaAttr[intersectFeaArea[max(intersectFeaArea.keys())]][1]) and temp_afiI > 0:
        countunderI.append(1)
    elif dictReferPolyAttr[row[indDIR]][1].contains(dictFeaAttr[intersectFeaArea[max(intersectFeaArea.keys())]][1]) and temp_afiI < 0:
        countoverI.append(1)
    #### the same code of the second logic ####
    
    if np.all(ReferPolyWeights == 1):
        pass
    else:
        rShapeArea = [dictReferPolyAttr[i][0] for i in dictReferPolyAttr if i in listfeaindDIR]
        rShapeArea.reverse()
        ReferPolyWeights = rShapeArea
    underMerging = np.average(underMergingI, axis = 0, weights = ReferPolyWeights); overMerging = np.average(overMergingI, axis = 0, weights = ReferPolyWeights)
    afi = np.average(afiI, weights = ReferPolyWeights);
    simsize = np.average(simsizeI, weights = ReferPolyWeights)
    rasub = np.average(rasubI, weights = ReferPolyWeights)
    rasuper = np.average(rasuperI, weights = ReferPolyWeights)
    qr = np.average(qrI, weights = ReferPolyWeights)
    oversegmentation = np.average(oversegmentationI, weights = ReferPolyWeights)
    undersegmentation = np.average(undersegmentationI, weights = ReferPolyWeights)
    rpsub = np.average(rpsubI, weights = ReferPolyWeights)
    rpsuper = np.average(rpsuperI, weights = ReferPolyWeights)
    dW = np.average(dIW, weights = ReferPolyWeights)
    d = np.average(dI, weights = ReferPolyWeights)
    mergesum = np.average(mergesumI, weights = ReferPolyWeights)
    MW = np.average(MIW, weights = ReferPolyWeights)
    M = np.average(MI, weights = ReferPolyWeights)
    zh1W = np.average(zh1IW, weights = ReferPolyWeights)
    zh1 = np.average(zh1I, weights = ReferPolyWeights)
    zh2W = np.average(zh2IW, weights = ReferPolyWeights)
    zh2 = np.average(zh2I, weights = ReferPolyWeights)
    return underMerging, afi, simsize, rasub, rasuper, qr, oversegmentation, undersegmentation, rpsub, rpsuper, d, dW, mergesum, M, MW, zh1, zh1W, zh2, zh2W
if __name__ == '__main__':
    dictReferPolyAttr, dictReferPointAttr, ReferPolyWeights = comborand_sample(sample_rate = 1)
    ReferPolyWeights = [1,2,3,4]
    c = combosegment_evaluation('1and1Scale100samwv2smallv6.shp', 0.5, dictReferPolyAttr, dictReferPointAttr, ReferPolyWeights)