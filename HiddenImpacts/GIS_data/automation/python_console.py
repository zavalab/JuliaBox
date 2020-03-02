#return all layers in the project:
layers = qgis.core.QgsProject.instance().layerTreeRoot().layerOrder()

#return all shown layers in the project:
iface.mapCanvas().layers()

#return selected layer in the project:
iface.activeLayer()

#set new selected layer:
iface.setActiveLayer(layer object)

#show/hide selected layer:
a = iface.actionShowSelectedLayers() -> a is a QAction object
a.trigger() -> to execute the action

#read project:
project_path = 'project.qgs'
template_path = 'temp.qpt'
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
document = QDomDocument()  from qgis.PyQt.QtXml import QDomDocument
document.setContent(template_content)

#open the layout:
composition = QgsLayout(this_proj)
composition.loadFromTemplate(document, QgsReadWriteContext())

#export as image:
dpi = 300;
exp = QgsLayoutExporter(composition)
exp.renderPageToImage(0,QSize(300, 300),dpi).save('out.bmp','bmp')

#accessing processing tools:
from processing.core.Processing import Processing
Processing.initialize()
from processing.tools import *
layerInput = QgsVectorLayer('test.shp','test','ogr')
general.run('qgis:explodelines', layerInput,'temp.shp')/processing.run(...)

check the help doc for algorithms:
processing.algorithmHelp("qgis:union")

#do interpolation
import qgis.analysis
layer = layers[0]
layer_data = QgsInterpolator.LayerData()
layer_data.source = layer
layer_data.interpolationAttribute = 3
idw_interpolator = QgsIDWInterpolator([layer_data])
idw_interpolator.setDistanceCoefficient(2)
export_path = 'test.asc'
rect = layer.extent()  or use layers[1].extent()
ncol = 500
nrow = 500
output = QgsGridFileWriter(idw_interpolator,export_path,rect,ncol, nrows)
output.writeFile()
iface.addRasterLayer(export_path, "some name")

#do clip:
params = {'INPUT': layers[0],
                  'MASK': layers[2],
                  'NODATA': 255,
                  'ALPHA_BAND': False,
                  'CROP_TO_CUTLINE': True,
                  'KEEP_RESOLUTION': False,
                  'OPTIONS': 'COMPRESS=LZW',
                  'DATA_TYPE': 6,
                  'OUTPUT': 'test.tif',
                  }
processing.run("gdal:cliprasterbymasklayer",params)
iface.addRasterLayer("test.tif", "some name")

#change raster properties:
from PyQt5.QtCore import *
from PyQt5.QtGui import *

layer = iface.activeLayer()

renderer = layer.renderer()
provider = layer.dataProvider()
extent = layer.extent()

ver = provider.hasStatistics(1, QgsRasterBandStats.All)

stats = provider.bandStatistics(1, QgsRasterBandStats.All,extent, 0)

if ver is not False:
    print("minimumValue = ", stats.minimumValue)
    print("maximumValue = ", stats.maximumValue)

minv= stats.minimumValue
maxv = stats.maximumValue
step = (maxv-minv)/4
colDic = {'red':'#d7191c', 'orange':'#fdae61', 'yellow':'#ffffbf', 'green':'#abdda4', 'blue':'#2b83ba'}
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

#move the layer to a place:
layer = layers[0]
root = QgsProject.instance().layerTreeRoot()
mylayer = root.findLayer(layer.id())
myClone = mylayer.clone()
parent = mylayer.parent()
parent.insertChildNode(-1, myClone)
parent.removeChildNode(mylayer)

#set markers color and size for point vector layers:
symbol = layer.renderer().symbols(QgsRenderContext())
symbol = symbol[0]
symbol.setColor(QColor('#000000'))
symbol.setSize(0.5)
qgis.utils.iface.mapCanvas().refreshAllLayers()

#set color and lind width for polygon vector layers:
layer = iface.activeLayer()
from PyQt5 import QtGui
myVectorLayer = layer
myRenderer  = layer.renderer()
if layer.wkbType() == 3 or layer.wkbType() == 6:
    mySymbol1 = QgsFillSymbol.createSimple({'color':'#ffffff',
                                            'color_border':'black',
                                            'width_border':'.26'})
    myRenderer.setSymbol(mySymbol1)

layer.triggerRepaint()
iface.mapCanvas().refreshAllLayers()
