from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1024, 668]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [127.5, 127.5, 127.5]
renderView1.UseToneMapping = 1
renderView1.UseAmbientOcclusion = 1
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-410.7422316343786, 127.5, 106.14242636025794]
renderView1.CameraFocalPoint = [127.5, 127.5, 127.5]
renderView1.CameraViewUp = [-8.796929262553182e-18, 1.0, 2.216955407797707e-16]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 220.83647796503186
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.Shadows = 1
renderView1.AmbientSamples = 1
renderView1.SamplesPerPixel = 4
renderView1.Backgroundmode = 'Backplate'
renderView1.EnvironmentalBG = [0.9975585564965286, 0.9975585564965286, 0.9975585564965286]
renderView1.OSPRayMaterialLibrary = materialLibrary1

layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(849, 554)

SetActiveView(renderView1)

# processing pipeline
data = XMLImageDataReader(registrationName='ctBones.vti', FileName=['data/ctBones.vti'])
data.PointArrayStatus = ['Scalars_']
data.TimeArray = 'None'

# create a new 'TTK TopologicalSimplificationByPersistence'
tTKTopologicalSimplificationByPersistence1 = TTKTopologicalSimplificationByPersistence(registrationName='TTKTopologicalSimplificationByPersistence1', Input=data)
tTKTopologicalSimplificationByPersistence1.InputArray = ['POINTS', 'Scalars_']
tTKTopologicalSimplificationByPersistence1.PersistenceThreshold = 155.0
tTKTopologicalSimplificationByPersistence1.ThresholdIsAbsolute = 1

# create a new 'TTK ExTreeM'
tTKExTreeM1 = TTKExTreeM(registrationName='TTKExTreeM1', Input=tTKTopologicalSimplificationByPersistence1)
tTKExTreeM1.ScalarArray = ['POINTS', 'Scalars_']

# find source
tTKExTreeM1_1 = FindSource('TTKExTreeM1')

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=OutputPort(tTKExTreeM1_1,2))
calculator3.Function = 'isSplitLeaf=1'

# create a new 'Iso Volume'
isoVolume2 = IsoVolume(registrationName='IsoVolume2', Input=calculator3)
isoVolume2.InputScalars = ['POINTS', 'Result']
isoVolume2.ThresholdRange = [0.5, 1.0]

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=isoVolume2)

# create a new 'Tetrahedralize'
tetrahedralize2 = Tetrahedralize(registrationName='Tetrahedralize2', Input=extractSurface1)

# create a new 'TTK ConnectedComponents'
tTKConnectedComponents1 = TTKConnectedComponents(registrationName='TTKConnectedComponents1', Segmentation=tetrahedralize2)
tTKConnectedComponents1.FeatureMask = ['POINTS', 'None']
tTKConnectedComponents1.SegmentationSize = 1

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=tTKConnectedComponents1)
threshold1.Scalars = ['POINTS', 'ComponentSize']
threshold1.LowerThreshold = 0.49999999999999994
threshold1.UpperThreshold = 50000.0
threshold1.ThresholdMethod = 'Above Upper Threshold'

# create a new 'TTK GeometrySmoother'
tTKGeometrySmoother1 = TTKGeometrySmoother(registrationName='TTKGeometrySmoother1', Input=threshold1)
tTKGeometrySmoother1.IterationNumber = 10
tTKGeometrySmoother1.InputMaskField = ['POINTS', 'Result']

# --- set up visualization ---

# show data from tTKTopologicalSimplificationByPersistence1
tTKTopologicalSimplificationByPersistence1Display = Show(tTKTopologicalSimplificationByPersistence1, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'Scalars_'
scalars_LUT = GetColorTransferFunction('Scalars_')
scalars_LUT.RGBPoints = [0.0, 1.0, 0.8, 0.64, 124.28154754638672, 0.773, 0.3, 0.3, 127.5, 0.865, 0.865, 0.865, 255.0, 0.23137254902, 0.298039215686, 0.752941176471]
scalars_LUT.ColorSpace = 'RGB'
scalars_LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Scalars_'
scalars_PWF = GetOpacityTransferFunction('Scalars_')
scalars_PWF.Points = [0.0, 0.0, 0.5, 0.0, 14.854368209838867, 0.0, 0.5, 0.0, 89.12621307373047, 0.2192513346672058, 0.5, 0.0, 127.25242614746094, 0.0, 0.5, 0.0, 255.0, 0.0, 0.5, 0.0]
scalars_PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
tTKTopologicalSimplificationByPersistence1Display.Representation = 'Volume'
tTKTopologicalSimplificationByPersistence1Display.ColorArrayName = ['POINTS', 'Scalars_']
tTKTopologicalSimplificationByPersistence1Display.LookupTable = scalars_LUT
tTKTopologicalSimplificationByPersistence1Display.SelectTCoordArray = 'None'
tTKTopologicalSimplificationByPersistence1Display.SelectNormalArray = 'None'
tTKTopologicalSimplificationByPersistence1Display.SelectTangentArray = 'None'
tTKTopologicalSimplificationByPersistence1Display.OSPRayScaleArray = 'Scalars_'
tTKTopologicalSimplificationByPersistence1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKTopologicalSimplificationByPersistence1Display.SelectOrientationVectors = 'None'
tTKTopologicalSimplificationByPersistence1Display.ScaleFactor = 25.5
tTKTopologicalSimplificationByPersistence1Display.SelectScaleArray = 'Scalars_'
tTKTopologicalSimplificationByPersistence1Display.GlyphType = 'Arrow'
tTKTopologicalSimplificationByPersistence1Display.GlyphTableIndexArray = 'Scalars_'
tTKTopologicalSimplificationByPersistence1Display.GaussianRadius = 1.2750000000000001
tTKTopologicalSimplificationByPersistence1Display.SetScaleArray = ['POINTS', 'Scalars_']
tTKTopologicalSimplificationByPersistence1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKTopologicalSimplificationByPersistence1Display.OpacityArray = ['POINTS', 'Scalars_']
tTKTopologicalSimplificationByPersistence1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKTopologicalSimplificationByPersistence1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKTopologicalSimplificationByPersistence1Display.PolarAxes = 'PolarAxesRepresentation'
tTKTopologicalSimplificationByPersistence1Display.ScalarOpacityUnitDistance = 1.7320508075688774
tTKTopologicalSimplificationByPersistence1Display.ScalarOpacityFunction = scalars_PWF
tTKTopologicalSimplificationByPersistence1Display.OpacityArrayName = ['POINTS', 'Scalars_']
tTKTopologicalSimplificationByPersistence1Display.Shade = 1
tTKTopologicalSimplificationByPersistence1Display.IsosurfaceValues = [127.5]
tTKTopologicalSimplificationByPersistence1Display.SliceFunction = 'Plane'
tTKTopologicalSimplificationByPersistence1Display.Slice = 127

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKTopologicalSimplificationByPersistence1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 255.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKTopologicalSimplificationByPersistence1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 255.0, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
tTKTopologicalSimplificationByPersistence1Display.SliceFunction.Origin = [127.5, 127.5, 127.5]

# show data from tTKGeometrySmoother1
tTKGeometrySmoother1Display = Show(tTKGeometrySmoother1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'ComponentId'
componentIdLUT = GetColorTransferFunction('ComponentId')
componentIdLUT.InterpretValuesAsCategories = 1
componentIdLUT.AnnotationsInitialized = 1
componentIdLUT.RGBPoints = [581236.0, 0.231373, 0.298039, 0.752941, 1247494.5, 0.865003, 0.865003, 0.865003, 1913753.0, 0.705882, 0.0156863, 0.14902]
componentIdLUT.NanColor = [0.8196078431372549, 0.8431372549019608, 0.8588235294117647]
componentIdLUT.ScalarRangeInitialized = 1.0
componentIdLUT.Annotations = ['581236', '581236', '1913379', '1913379', '1913381', '1913381', '1913515', '1913515', '1913525', '1913525', '1913598', '1913598', '1913653', '1913653', '1913705', '1913705', '1913712', '1913712', '1913744', '1913744', '1913753', '1913753']
componentIdLUT.ActiveAnnotatedValues = ['581236', '1913379', '1913381', '1913515', '1913525', '1913598', '1913653', '1913705', '1913712', '1913744', '1913753']
componentIdLUT.IndexedColors = [0.8901960784313725, 0.8549019607843137, 0.788235294117647, 0.996078431372549, 0.8509803921568627, 0.5568627450980392, 1.0, 1.0, 0.8313725490196079, 0.996078431372549, 0.8509803921568627, 0.5568627450980392, 1.0, 1.0, 0.8313725490196079, 1.0, 1.0, 0.8313725490196079, 0.996078431372549, 0.8509803921568627, 0.5568627450980392, 0.996078431372549, 0.8509803921568627, 0.5568627450980392, 1.0, 1.0, 0.8313725490196079, 1.0, 0.9686274509803922, 0.7372549019607844, 0.996078431372549, 0.8509803921568627, 0.5568627450980392]
componentIdLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# get opacity transfer function/opacity map for 'ComponentId'
componentIdPWF = GetOpacityTransferFunction('ComponentId')
componentIdPWF.Points = [581236.0, 0.0, 0.5, 0.0, 1913753.0, 1.0, 0.5, 0.0]
componentIdPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
tTKGeometrySmoother1Display.Representation = 'Surface'
tTKGeometrySmoother1Display.ColorArrayName = ['POINTS', 'ComponentId']
tTKGeometrySmoother1Display.LookupTable = componentIdLUT
tTKGeometrySmoother1Display.SelectTCoordArray = 'None'
tTKGeometrySmoother1Display.SelectNormalArray = 'None'
tTKGeometrySmoother1Display.SelectTangentArray = 'None'
tTKGeometrySmoother1Display.OSPRayScaleArray = 'Result'
tTKGeometrySmoother1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKGeometrySmoother1Display.SelectOrientationVectors = 'None'
tTKGeometrySmoother1Display.ScaleFactor = 23.60305922329426
tTKGeometrySmoother1Display.SelectScaleArray = 'Result'
tTKGeometrySmoother1Display.GlyphType = 'Arrow'
tTKGeometrySmoother1Display.GlyphTableIndexArray = 'Result'
tTKGeometrySmoother1Display.GaussianRadius = 1.180152961164713
tTKGeometrySmoother1Display.SetScaleArray = ['POINTS', 'Result']
tTKGeometrySmoother1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKGeometrySmoother1Display.OpacityArray = ['POINTS', 'Result']
tTKGeometrySmoother1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKGeometrySmoother1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKGeometrySmoother1Display.PolarAxes = 'PolarAxesRepresentation'
tTKGeometrySmoother1Display.ScalarOpacityFunction = componentIdPWF
tTKGeometrySmoother1Display.ScalarOpacityUnitDistance = 2.550379689216199
tTKGeometrySmoother1Display.OpacityArrayName = ['POINTS', 'Result']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKGeometrySmoother1Display.ScaleTransferFunction.Points = [0.49999999999999994, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKGeometrySmoother1Display.OpacityTransferFunction.Points = [0.49999999999999994, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

Render()

#save screenshot
WriteImage("output/test.png")

