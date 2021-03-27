# Created by: Ben Gannon
# Modified by: Stephanie Mueller
# Created on: 09/16/2016
# Last Updated: 01/21/2021

'''
This script is designed to create six products related to treatment
feasibility and cost. The three fuel treatment options currently being
considered are mechanical only, Rx fire only, and complete (mechanical
and prescribed (Rx) fire).

Each treatment type outputs a cost and feasibility(binary) raster.

# Thin only
# 1) Cost as function of distance from roads and slope
# 2) Feasibility (restricted from wilderness and Upper Tier roadless)
#
# Rx fire only
# 3) Costs assumed to be uniform based on distance from WUI.
#    Based on communication with local fuels and fire planners (Bryan Karchut (ARP) and
#    James White (CLRD)), current cost is close to $1000/ac when > 250 m
#    from structures. Also includes option to include area within 250 m
#    WUI buffer at an increased cost of $2000/ac
# 4) Feasibility (restricted from within 250 m of structures, and frequent fire EVT types)
# 
# Thin and Rx fire
# 5) Cost assumed to be the sum of thin only and Rx fire only costs
# 6) Feasibility assumed to be the area of thin only
'''

#-> Import modules
import arcpy, os, sys, traceback, math, shutil, csv
from arcpy.sa import *

start = datetime.datetime.now()

#-> Check out spatial analyst extension
arcpy.CheckOutExtension('spatial')
arcpy.env.overwriteOutput = True

#-> function(s)
def checkPath(path):
    try: os.makedirs(path)
    except: pass

#-> Create workspace
plist = sys.argv[0].split("\\")                         # Grab script folder pathname
mpath = "\\".join(plist[0:(len(plist)-2)])
arcpy.env.workspace = mpath + r'\scripts' 
arcpy.env.scratchWorkspace = mpath + r'\scripts\INTERMEDIATE'

#--> Create variables

# Input rasters
slope = r'\INPUT\SPATIAL\RAW_LANDFIRE\us_slp2010'              
cc = r'\INPUT\SPATIAL\RAW_LANDFIRE\cc'
evt = r'\INPUT\SPATIAL\RAW_LANDFIRE\evt'

# Input feature classes
extent = r'INPUT\SPATIAL\VECTOR_INPUT.gdb\raster_extent'
nothin = r'\INPUT\SPATIAL\VECTOR_INPUT.gdb\exclude_thinning'    
roads = r'\INPUT\SPATIAL\VECTOR_INPUT.gdb\Trans_RoadSegment'
WUI_structures = r'\INPUT\SPATIAL\VECTOR_INPUT.gdb\MS_Bing_structures'

# Input reclass table
Rx_evt_feas = r'\INPUT\Rx_evt_feasibility.csv'             

# Output
constraints_dir = arcpy.env.workspace + r'\INPUT\SPATIAL\Constraints'
checkPath(constraints_dir)

mocost_out = constraints_dir + r'\mocost.tif'
mofeas_out = constraints_dir + r'\mofeas.tif'
Rxocost_out = constraints_dir + r'\Rxocost.tif'
Rxofeas_out = constraints_dir + r'\Rxofeas.tif'
cfeas_out = constraints_dir + r'\cfeas.tif'
ccost_out = constraints_dir + r'\ccost.tif'

# Model variables
MaxCost = 10000 # maximum mechanical cost $/ac                                   
BaseCost = 2500 # base mechanical cost $/ac
RxFireCost = 1000 # $/ac - Bryan Karchut (ARP personal communication)
IncludeRXWUI = 1 # 1 = include (feasible adjacent to WUI), 0 = not included (not feasible within < 250 m from WUI)
RxFireWUICost = 3000 # Rx Fire cost in WUI (< 250 m from WUI) $/ac
slopeTH = 40 # percent - math.atan(0.4)*(180/math.pi) 
slopeMAX = 200 # percent  
distTH = 800 # meters            
distMAX = 6400 # meters
ccTH = 10 # Canopy cover of <10% (To limit mechanical treatments to forested fuels)

#######-------> Main body of analysis
print('###---> Modeling treatment costs and feasibilities')

# Set snap raster and extent
arcpy.env.snapRaster = cc
arcpy.env.extent = cc

# Convert extent to raster to avoid file creations by zonal stats and tabulate area
try:
    # Add zone field and set to 1
    arcpy.AddField_management(extent,'ZONE','SHORT')
    arcpy.CalculateField_management(extent,'ZONE','1','PYTHON_9.3')
    # Convert to raster
    extent_r = arcpy.env.scratchWorkspace + r'\extent_r'
    arcpy.FeatureToRaster_conversion(extent,'ZONE',extent_r,30)
    print('Converted extent to raster for zonal stats and tabulate area')
except Exception as err:
    print(err.args[0])

print

#####-----> Mechanical treatment costs
print('###---> Working on mechanical only treatment costs')

###---> Deal with costs associated with distance from roads

# Convert roads to raster and then to binary raster
try:
    roads_r = arcpy.env.scratchWorkspace + r'\roads' 
    arcpy.PolylineToRaster_conversion(roads,'OBJECTID',roads_r,'MAXIMUM_LENGTH','',30)
    rbin = Con(Raster(roads_r) > 0,1,0)
    print('Converted roads feature class to raster and converted to binary')
except Exception as err:
    print(err.args[0])

# Calculate Euclidean distance from roads binary raster
try:
    Euc_Dist = EucDistance(rbin)
    print('Calculated Euclidean distance')
except Exception as err:
    print(err.args[0])

# Calculate road-distance additional costs
# Note assumed maximum distance of 4 miles (6400 m)
try:
    m1 =(MaxCost-BaseCost)/(distMAX-distTH)
    rcost = Con(Euc_Dist <= distTH,0,(Euc_Dist-distTH)*m1)
    print('Calculated road-distance additional cost')
except Exception as err:
    print(err.args[0])

###---> Deal with costs associated with slope (NOTE: slope input is in deg)

# Calculate slope additional costs
try:
    # Convert slope to percent
    test = arcpy.env.scratchWorkspace + r'\PerSlope'
    PerSlope = Tan(Float(Raster(slope))*(math.pi/180))*100
    PerSlope.save(test)
    # Calculate costs
    m2 = (MaxCost-BaseCost)/(slopeMAX-slopeTH)
    scost = Con(PerSlope <= slopeTH,0,(PerSlope-slopeTH)*m2)
    print('Calculated slope additional cost')
except Exception as err:
    print(err.args[0])
   
###---> Combine base, road, and slope costs
try:
    mocost = rcost + scost + BaseCost
    # Limit to max cost
    mocost = Con(mocost > MaxCost,MaxCost,mocost)
    mocost.save(mocost_out)
    print('Combined base, road, and slope costs')
except Exception as err:
    print(err.args[0])

# Print final raster statistics to shell
MCmin = arcpy.GetRasterProperties_management(mocost,'MINIMUM')
MCmean = arcpy.GetRasterProperties_management(mocost,'MEAN')
MCmax = arcpy.GetRasterProperties_management(mocost,'MAXIMUM')

print
print('Final Raster Attributes:')
print('Minimum raster value = ' + str(MCmin.getOutput(0)))
print('Mean raster value = ' + str(MCmean.getOutput(0)))
print('Maximum raster value = ' + str(MCmax.getOutput(0)))
print

#####-----> Mechanical treatment feasibility
print('###---> Working on mechanical only treatment feasibility')
# Merge feature layers
try:
    # Exlude Thinning
    nothin_cl = arcpy.env.scratchWorkspace + r'\nothin.shp'
    arcpy.Clip_analysis(nothin,extent,nothin_cl)

except Exception as err:
    print(err.args[0])

# Convert to binary raster
try:
    # Assign all polygons the NOTREAT value
    nothin_fl = 'nothin_fl'
    arcpy.MakeFeatureLayer_management(nothin_cl,nothin_fl)
    arcpy.AddField_management(nothin_fl,'NOTREAT','SINGLE')
    arcpy.CalculateField_management(nothin_fl,'NOTREAT',0,'PYTHON_9.3')
    # Prepare LandFire extent [to have background value]
    extent_fl = 'lf_extent_fl'
    arcpy.MakeFeatureLayer_management(extent,extent_fl)    
    arcpy.AddField_management(extent_fl,'NOTREAT','SINGLE')
    arcpy.CalculateField_management(extent_fl,'NOTREAT',1,'PYTHON_9.3')
    # Update with LandFire extent
    update_fc = arcpy.env.scratchWorkspace + r'\update_fc.shp'
    arcpy.Update_analysis(extent_fl,nothin_fl,update_fc)
    # Convert to raster
    mgmtrest = arcpy.env.scratchWorkspace + r'\mgmtrest'
    arcpy.PolygonToRaster_conversion(update_fc,'NOTREAT',mgmtrest,'','',30)
    print('Made land management restrictions raster')
except Exception as err:
    print(err.args[0])

# Apply canopy cover threshold
try:
    ccBin = Con(Raster(cc) >= ccTH,1,0)
    print('Created canopy cover threshold raster')
except Exception as err:
    print(err.args[0])
    
###---> Assemble sub-factors
try:
    mofeas = ccBin*mgmtrest
    mofeas.save(mofeas_out)
    print('Created and saved mechanical only feasibility raster')
except Exception as err:
    print(err.args[0])

# Summarize results and print to shell
try:
    # Summarize using tabulate area
    temp_t = arcpy.env.scratchWorkspace + r'\temp_t'
    if arcpy.Exists(temp_t):
        arcpy.Delete_management(temp_t)
    TabulateArea(extent_r,'VALUE',mofeas,'VALUE',temp_t)
    # Access results with search cursor
    fields_S = ['VALUE_0','VALUE_1']    # Fields for search cursor
    with arcpy.da.SearchCursor(temp_t,fields_S) as cursor:
        for row in cursor:
            untreatable = row[0]
            treatable = row[1]
            total = row[0] + row[1]
    del cursor
    print
    print('Final Raster Attributes:')
    print('Percent treatable:   ' + str((treatable/total)*100))
    print('Percent untreatable: ' + str((untreatable/total)*100))
    print
except Exception as err:
    print(err.args[0])

#####-----> Prescribed fire only treatment costs
print('###---> Working on prescribed fire only treatment costs')

###---> Do interpolation with Inverse Distance Weighting (IDW)
# Convert structures point data to raster
try:
    scount = arcpy.env.scratchWorkspace + r'\scount'
    arcpy.PointToRaster_conversion(WUI_structures,'',scount,'COUNT','',30)
    print('Generated raster of structure counts')
except Exception as err:
    print(err.args[0])

# Calculate Euclidean Distance
try:
    dist = EucDistance(Raster(scount))
    print('Calculated Euclidean Distance')
except Exception as err:
    print(err.args[0])

# Calculate costs based on distance to WUI
try:
    Rxocost = Con(dist > 1000,RxFireCost,RxFireWUICost)
    Rxocost.save(Rxocost_out)
    print('Created prescribed fire cost raster')
except Exception as err:
    print(err.args[0])

# Calculate WUI Rx additional cost

# Print final raster statistics to shell
Rxmin = arcpy.GetRasterProperties_management(Rxocost,'MINIMUM')
Rxmean = arcpy.GetRasterProperties_management(Rxocost,'MEAN')
Rxmax = arcpy.GetRasterProperties_management(Rxocost,'MAXIMUM')

print
print('Final Raster Attributes:')
print('Minimum raster value = ' + str(Rxmin.getOutput(0)))
print('Mean raster value = ' + str(Rxmean.getOutput(0)))
print('Maximum raster value = ' + str(Rxmax.getOutput(0)))
print
 
#####-----> Prescribed fire treatment feasibility
print('###---> Working on Rx fire only treatment feasibility')

# Classify into Rx fire feasibility raster
try:
    if IncludeRXWUI == 0:
        srest = Con(dist > 250,1,0)
        print('Created structure restrictions raster > 250 m from WUI')
    else:
        srest = Con(dist > 0,1,0)
        print('Created raster without WUI restrictions')
except Exception as err:
    print(err.args[0])

###---> Execute Reclassify of EVT
try:
    Rxeco = ReclassByTable(evt, Rx_evt_feas, "VALUE", "VALUE", "FEASIBLE", "NODATA")  
    Rxeco.save(arcpy.env.scratchWorkspace + "/Rxeco")
    print('Remapped EVT to ecologically-appropriate Rx application categories')
except Exception as err:
    print(err.args[0])

###---> Assemble sub-factors
try:
    Rxofeas = srest*Rxeco
    Rxofeas.save(Rxofeas_out)
    print('Created and saved prescribed fire feasibility raster')
except Exception as err:
    print(err.args[0])

# Summarize results and print to shell
try:
    # Summarize using tabulate area
    temp_t = arcpy.env.scratchWorkspace + r'\temp_t'
    if arcpy.Exists(temp_t):
        arcpy.Delete_management(temp_t)
    TabulateArea(extent_r,'VALUE',Rxofeas,'VALUE',temp_t)
    # Access results with search cursor
    fields_S = ['VALUE_0','VALUE_1']    # Fields for search cursor
    with arcpy.da.SearchCursor(temp_t,fields_S) as cursor:
        for row in cursor:
            untreatable = row[0]
            treatable = row[1]
            total = row[0] + row[1]
    del cursor
    print
    print('Final Raster Attributes:')
    print('Percent treatable:   ' + str((treatable/total)*100))
    print('Percent untreatable: ' + str((untreatable/total)*100))
    print
except Exception as err:
    print(err.args[0])

#####-----> Mechanical and Rx fire treatment costs
print('###---> Working on mechanical and Rx fire treatment costs')

# Generate constant raster as place holder
try:
    ccost = mocost + Rxocost
    ccost.save(ccost_out)
    print('Created mechanical and Rx fire cost raster')
except Exception as err:
    print(err.args[0])

# Print final raster statistics to shell
mRxmin = arcpy.GetRasterProperties_management(ccost,'MINIMUM')
mRxmean = arcpy.GetRasterProperties_management(ccost,'MEAN')
mRxmax = arcpy.GetRasterProperties_management(ccost,'MAXIMUM')

print
print('Final Raster Attributes:')
print('Minimum raster value = ' + str(mRxmin.getOutput(0)))
print('Mean raster value = ' + str(mRxmean.getOutput(0)))
print('Maximum raster value = ' + str(mRxmax.getOutput(0)))
print

#####-----> Combined treatment feasibility
print('###---> Working on combined treatment feasibility')

# Combine mechanical and prescribed fire feasibility rasters
try:
    #cmRxfeas = mofeas
    cfeas = mofeas
    cfeas.save(cfeas_out)
    print('Created combined feasibility raster')
except Exception as err:
    print(err.args[0])

# Summarize results and print to shell
try:
    # Summarize using tabulate area
    temp_t = arcpy.env.scratchWorkspace + r'\temp_t'
    if arcpy.Exists(temp_t):
        arcpy.Delete_management(temp_t)
    TabulateArea(extent_r,'VALUE',cfeas,'VALUE',temp_t)
    # Access results with search cursor
    fields_S = ['VALUE_0','VALUE_1']    # Fields for search cursor
    with arcpy.da.SearchCursor(temp_t,fields_S) as cursor:
        for row in cursor:
            untreatable = row[0]
            treatable = row[1]
            total = row[0] + row[1]
    del cursor
    print
    print('Final Raster Attributes:')
    print('Percent treatable:   ' + str((treatable/total)*100))
    print('Percent untreatable: ' + str((untreatable/total)*100))
    print
except Exception as err:
    print(err.args[0])

###---> Clean up the scratch workspace

# Set the workspace to the SCRATCH folder so you can grab the rasters
arcpy.env.workspace = arcpy.env.scratchWorkspace
rasters = arcpy.ListRasters('*','ALL')
features = arcpy.ListFeatureClasses('*','ALL')
tables = arcpy.ListTables('*','ALL')

# Delete all the rasters, feature classes, and tables
try:
    for raster in rasters:
        arcpy.Delete_management(raster)
    for feature in features:
        arcpy.Delete_management(feature)
    for table in tables:
        arcpy.Delete_management(table)
except Exception as err:
    print(err.args[0])

print "Script complete! Run Time: %s\n\n" %(datetime.datetime.now() - start)

