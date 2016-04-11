# -*- coding: cp936 -*-
def export_fea_names(curDir):
    import arcpy
    feature_names = []     # store the input feature names
    workspace = curDir + '\\results'
    walk = arcpy.da.Walk(workspace, datatype="FeatureClass", type="Polygon")
    # set the corresponding variables to store the corresponding files 
    for dirpath, dirnames, filenames in walk:
        for filename in filenames:
            if filename.find('.'):
                feature_names.append(filename)
    return feature_names

    
if __name__ == '__main__':
    import os
    feature_names = export_fea_names(os.getcwd())