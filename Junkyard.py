# This is a junkyard for old methods or code that I want to keep around just in case.

def CalcChannelWidthAdvancedLic(nodexy, rb, lb):
    rb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, rb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(rb_near_result, ['NEAR_DIST']) as rows:
            row = rows.next()
            rb_distance = row[0]    
    
    lb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, lb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(lb_near_result, ['NEAR_DIST']) as rows:
                row = rows.next()
                lb_distance = row[0]    
    return(rb_distance + lb_distance)

def GetTopoAngles2(NodeDict, streamID, EleRaster, azimuths, azimuthdisdict, con_z_to_m):
    """This gets the maximum shade angle for each stream and outputs the data as a list. 
    This version generates an array for each stream. It is faster but may crash if the arrays are large"""
    
    TopoList = []
    CoordList = []
    
    cellsizeX = arcpy.Describe(EleRaster).meanCellWidth
    cellsizeY = arcpy.Describe(EleRaster).meanCellHeight
    
    # calculate the buffer distance (in EleRaster spatial units) to add to the raster bounding box when extracting to an array
    buffer = cellsizeX * 2    
    
    Nodes = NodeDict.keys()
    Nodes.sort()
    
    for nodeID in Nodes:
        origin_x = NodeDict[nodeID]["POINT_X"]
        origin_y = NodeDict[nodeID]["POINT_Y"]      
        for a in azimuths:
            for searchdistance in azimuthdisdict[a]:
                # Calculate the x and y coordinate of the landcover sample location in the spatial units of the inputFC                
                topoX = (searchdistance * con_from_m * sin(radians(a))) + origin_x
                topoY = (searchdistance * con_from_m * cos(radians(a))) + origin_y

                # Add the all the coordinates to the list
                CoordList.append([topoX, topoY, nodeID, a])
            
    # calculate lower left corner coordinate and nrows/cols for the bounding box
    # first transpose the list so x and y coordinates are in the same vector
    tlist = map(lambda *i: list(i), *CoordList)
      
    Xmin = min(tlist[0]) - buffer
    Ymin = min(tlist[1]) - buffer
    Ymax = max(tlist[1]) + buffer            
    ncols = (max(tlist[0]) + buffer - Xmin) / cellsizeX
    nrows = (Ymax - Ymin) / cellsizeY
    bbox_lower_left = arcpy.Point(Xmin, Ymin) # must be in raster map units
    bbox_upper_left = [Xmin, Ymax, cellsizeX, cellsizeY]
    nodata_to_value = -9999 / con_z_to_m
    
    del(tlist)
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    Zarry = arcpy.RasterToNumPyArray(EleRaster, bbox_lower_left, ncols, nrows, nodata_to_value)
    
    # convert array values to meters if needed
    Zarry = Zarry * con_z_to_m             
    
    # iterate through each node and azimuth
    for nodeID in NodeDict:
        for a in azimuths:
            CoordList2 = [row for row in CoordList if row[2] == nodeID and row[3] == a]
            TopoZList = []
            #xy = []
    
            for i in range(0,len(CoordList2)):
                xy = CoordToArray(CoordList2[i][0], CoordList2[i][1], bbox_upper_left)
                TopoZList.append(Zarry[xy[1], xy[0]])
                #xy.append(CoordToArray(CoordList2[i][0], CoordList2[i][1], bbox_upper_left))
                #TopoZList.append(1)
            
            TopoZarry = np.array(TopoZList)
            Disarry = azimuthdisdict[a] * con_to_m
            Shadearry = np.degrees(np.arctan((TopoZarry - TopoZarry[0]) / Disarry))
            # Take out the off raster samples
            naindex = np.where(TopoZarry < -9998)
            for x in naindex[0]: Shadearry[x] = -9999            
            # Find the max shade angle
            ShadeAngle = Shadearry.max()
            # array index at the max shade angle 
            arryindex = np.where(Shadearry==ShadeAngle)[0][0]
            ShadeZ = TopoZarry[arryindex]
            ZChange = ShadeZ - TopoZarry[0]
            ShadeDistance = azimuthdisdict[a][arryindex]
            SearchDistance = azimuthdisdict[a].max()
            ShadeAngle_X = (ShadeDistance * con_from_m * sin(radians(a))) + CoordList2[0][0]
            ShadeAngle_Y = (ShadeDistance * con_from_m * cos(radians(a))) + CoordList2[0][1]
            offRasterSamples = (TopoZarry > -9998).sum()
        
            TopoList.append([ShadeAngle_X, ShadeAngle_Y, ShadeAngle_X, ShadeAngle_Y, streamID, nodeID, a, ShadeAngle, ShadeZ, TopoZarry[0], ZChange, ShadeDistance, SearchDistance, offRasterSamples])

    return TopoList

def Downloadcomtypes():
    import urllib
    
    print("Downloading comtypes")
    url = "https://pypi.python.org/packages/source/c/comtypes/comtypes-1.1.1.zip"
    file_name = url.split('/')[-1]
    out_filepath = "C:/" + file_name
    urllib.urlretrieve (url, r"C:\comtypes-1.1.1.zip")

def GetLibPath():
    """Return location of ArcGIS type libraries as string.
    Based on work from Matt Wilkie and others
    https://bitbucket.org/maphew/canvec/src/eaf2678de06f/Canvec/Scripts/parco.py"""
    import _winreg

    # this looks for a file which should always be present in the ESRI  com library folder
    # first in the regkey (and directory) expected for v9x, then for v10
    # this method of testing is ugly, reeks of code smell, but works
    # after 2 hours of trying to make it happen elegantly I give up, for now. (mhw)

    keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\ArcGIS")
    esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"

    if not os.path.isfile(esriComDir+'\\esriSystemUI.olb'):
        keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\Desktop10.1")
        esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"
    
    if not os.path.isfile(esriComDir+'\\esriSystemUI.olb'):
            keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\Desktop10.2")
            esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"    

    return esriComDir
    
def GetModule(sModuleName):
    """Import ArcGIS module"""
    from comtypes.client import GetModule
    sLibPath = GetLibPath()
    GetModule(sLibPath + sModuleName)

def GetStandaloneModules():
    """Import commonly used ArcGIS libraries for standalone scripts"""
    GetModule("esriSystem.olb")
    GetModule("esriGeometry.olb")
    GetModule("esriCarto.olb")
    GetModule("esriDisplay.olb")
    GetModule("esriGeoDatabase.olb")
    GetModule("esriDataSourcesGDB.olb")
    GetModule("esriDataSourcesFile.olb")
    GetModule("esriOutput.olb")

def DistanceBetweenPoints(Xa,Ya,Xb,Yb):
    """Determines the distance between two points using the pythagorean theorem and returns distance between them in map units"""
    from math import sqrt, pow
    dist = math.sqrt(math.pow((Xa - Xb),2) + math.pow((Ya - Yb),2))
    return dist

def CalcAspect(dx, dy):
    
    RadAngle = atan2(dY,dX)
    aspect = (RadAngle * 180) / pi
    aspect2 = degrees(atan2(dY,dX))