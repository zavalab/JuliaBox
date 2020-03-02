# To do interpolation and generate maps using QGIS 3.2.2
# Require a QGIS project file (.qgz/.qhs), a base map (.shp),
# the study region boundary (.shp), node file with attributes (.csv),
# and a template to compose on (.qpt)
# Yicheng Hu 2018-Nov

from qgis.core import *
from qgis.gui import *
from qgis.analysis import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from qgis.PyQt.QtXml import QDomDocument
from processing.core.Processing import Processing
Processing.initialize()
from processing.tools import *

global_path = '/Users/asampat/research/BioGas/Code/Apoorva/hidden_impacts/new_data/GIS data/automation'
project_path = 'test.qgz' # .qgs also works
template_path = 'temp.qpt'
basemap_path = 'Dane_Municipalities_2015.shp'
boundary_path = 'AOI_Boundary_new.shp'  # the crs should coordinate with the csv
csvfile = '45' # name of the csv file without .csv

#open a project
this_proj = QgsProject.instance()
this_proj.read(project_path)
canvas = QgsMapCanvas()
bridge = QgsLayerTreeMapCanvasBridge(
    QgsProject.instance().layerTreeRoot(), canvas)
bridge.setCanvasLayers()

#read template:
template_file = open(template_path)
template_content = template_file.read()
template_file.close()
document = QDomDocument()
document.setContent(template_content)

#add base map and set up the properties
iface.addVectorLayer(basemap_path, 'base map', 'ogr');
base_map = iface.activeLayer()
myVectorLayer = base_map
myRenderer  = base_map.renderer()
if base_map.wkbType() == 3 or base_map.wkbType() == 6:  # for polygon vector layer
    mySymbol1 = QgsFillSymbol.createSimple({'color':'255,0,0,0', # transparent
                                            'color_border':'black', # stroke color
                                            'width_border':'.26'}) # stroke width
    myRenderer.setSymbol(mySymbol1)

base_map.triggerRepaint()
iface.mapCanvas().refreshAllLayers()

#add boundary map and set up the properties
iface.addVectorLayer(boundary_path, 'boundary', 'ogr');
boundary = iface.activeLayer()
from PyQt5 import QtGui
myVectorLayer = boundary
myRenderer  = boundary.renderer()
if boundary.wkbType() == 3 or boundary.wkbType() == 6:
    mySymbol1 = QgsFillSymbol.createSimple({'color':'255,0,0,0',
                                            'color_border':'black',
                                            'width_border':'.66'})
    myRenderer.setSymbol(mySymbol1)

boundary.triggerRepaint()
iface.mapCanvas().refreshAllLayers()

# add points and set up properties
uri='file://' + global_path +  csvfile + '.csv?crs=epsg:4326&delimiter=,&yField=Latitude&xField=Longitude'
layer = QgsVectorLayer(uri, csvfile, 'delimitedtext')
this_proj.addMapLayer(layer)
symbol = layer.renderer().symbols(QgsRenderContext())
symbol = symbol[0]
symbol.setColor(QColor('#000000')) # node dot color black
symbol.setSize(0.5) # node dot size
qgis.utils.iface.mapCanvas().refreshAllLayers()

# do interpolation
layer = iface.activeLayer()
layer_data = QgsInterpolator.LayerData()
layer_data.source = layer
layer_data.interpolationAttribute = 3 # the number of columns of the attribute
idw_interpolator = QgsIDWInterpolator([layer_data])
idw_interpolator.setDistanceCoefficient(2) # set distance coefficient here, 2 is default
prod = layer.fields()[layer_data.interpolationAttribute].name()
export_path = csvfile + '-' + prod + '.asc'
rect = boundary.extent()
ncol = 500 # number of columns
nrow = 500 # number of rows
output = QgsGridFileWriter(idw_interpolator,export_path,rect,ncol, nrow)
output.writeFile()
iface.addRasterLayer(export_path, csvfile + '-' + prod)
raster = iface.activeLayer()

# do clip
output = csvfile + '-' + prod + '.tif'
params = {'INPUT': raster,
                  'MASK': boundary,
                  'NODATA': 255,
                  'ALPHA_BAND': False,
                  'CROP_TO_CUTLINE': True,
                  'KEEP_RESOLUTION': False,
                  'OPTIONS': 'COMPRESS=LZW',
                  'DATA_TYPE': 6, # float number
                  'OUTPUT': output,
                  }
processing.run("gdal:cliprasterbymasklayer",params)
iface.addRasterLayer(output, csvfile + '-' + prod  + '-clipped')
this_proj.removeMapLayer(raster.id())
iface.mapCanvas().refreshAllLayers()

# set up ramp colors
layer = iface.activeLayer()

renderer = layer.renderer()
provider = layer.dataProvider()
extent = layer.extent()

ver = provider.hasStatistics(1, QgsRasterBandStats.All)

stats = provider.bandStatistics(1, QgsRasterBandStats.All,extent, 0)

if ver is not False:
    print("minimumValue = ", stats.minimumValue)
    print("maximumValue = ", stats.maximumValue)

minv= stats.minimumValue # change the price range here
maxv = stats.maximumValue
step = (maxv-minv)/4
colDic = {'red':'#d7191c', 'orange':'#fdae61', 'yellow':'#ffffbf', 'green':'#abdda4', 'blue':'#2b83ba'} # the ramp color is set up here: 5 colors
valueList =[minv, minv+step, minv+2*step, minv+3*step, maxv]
lst = [ QgsColorRampShader.ColorRampItem(valueList[0], QColor(colDic['red']), str(valueList[0])),
        QgsColorRampShader.ColorRampItem(valueList[1], QColor(colDic['orange']), str(valueList[1])),
        QgsColorRampShader.ColorRampItem(valueList[2], QColor(colDic['yellow']), str(valueList[2])),
        QgsColorRampShader.ColorRampItem(valueList[3], QColor(colDic['green']), str(valueList[3])),
        QgsColorRampShader.ColorRampItem(valueList[4], QColor(colDic['blue']), str(valueList[4]))]

myRasterShader = QgsRasterShader()
myColorRamp = QgsColorRampShader()
myColorRamp.setColorRampItemList(lst)
myColorRamp.setColorRampType(QgsColorRampShader.Interpolated)
myRasterShader.setRasterShaderFunction(myColorRamp)
myPseudoRenderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(),
                                                    layer.type(),
                                                    myRasterShader)
layer.setRenderer(myPseudoRenderer)
layer.triggerRepaint()

# move layer: move the clipped raster layer to the end
root = QgsProject.instance().layerTreeRoot()
mylayer = root.findLayer(layer.id())
myClone = mylayer.clone()
parent = mylayer.parent()
parent.insertChildNode(-1, myClone)
parent.removeChildNode(mylayer)

# export the map
#open the layout:
composition = QgsLayout(this_proj)
composition.loadFromTemplate(document, QgsReadWriteContext())

#export as image:
dpi = 300 # the resolution
exp = QgsLayoutExporter(composition)
exp.renderPageToImage(0,QSize(300, 300),dpi).save(csvfile + '-' + prod + '.bmp','bmp')
