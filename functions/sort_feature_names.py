import os 
import sys

def strfeacmp(featurename1, featurename2):
    import re
    pattern = re.compile(r'\d+')
    fea_nums1 = pattern.findall(featurename1)
    fea_nums2 = pattern.findall(featurename2)
    fea_nums1 = [int(i) for i in fea_nums1]
    fea_nums2 = [int(i) for i in fea_nums2]
    for i in range(0,len(fea_nums1)-2):
        if fea_nums1[i] == fea_nums2[i]:
            continue
        else:
            if fea_nums1[i] > fea_nums2[i]:
                return True
            else:
                return False
# this is insert sort
def FeatureQuickSort(arr):  
    for i in range(1,len(arr)):  
        j=i  
        while j>0 and strfeacmp(arr[j-1],arr[i]):  
            j-=1  
        arr.insert(j,arr[i])  
        arr.pop(i+1)  
        
if __name__ == '__main__':
    import arcpy
    curDir = os.getcwd()
    arcpy.env.workspace = curDir + '\\results'
    feature_names = []    
    workspace = curDir + '\\results'
    walk = arcpy.da.Walk(workspace, datatype="FeatureClass", type="Polygon")
    
    for dirpath, dirnames, filenames in walk:
        for filename in filenames:
            if filename.find('.'):
                feature_names.append(filename)
    FeatureQuickSort(feature_names)
    flag = strfeacmp('1and1Scale1samwv2PCAv3','1and1Scale6samwv2PCAv3')
    print repr(flag)
    