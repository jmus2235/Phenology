//EVI

// Add 30m DEM
var image = ee.Image('JAXA/ALOS/AW3D30_V1_1')
print(dsm)

// Select average DEM class
var elev = dsm.select('AVE')

// Select elevations greater than 500m
var elevGT500 = elev.gt(500)
elevGT500 = elevGT500.updateMask(elevGT500)
var mask = elev.gt(500)
//print(elevGT500)

// Add mask as a layer to display
//Map.addLayer(elevGT500, {palette: ['yellow']}, 'elevGT500')

//load the image collection
var L8EVI = ee.ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI') //Landsat 8 Collection 1 Tier 1 8-Day EVI Composite
    .filter(ee.Filter.dayOfYear(0, 365))
    .filterDate('2015-01-01', '2018-12-31')
    .filterBounds(BART);
print(L8EVI)

var L8TOA = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
    .filter(ee.Filter.dayOfYear(0, 365))
    .filterDate('2015-08-19', '2015-08-28')
    .filterBounds(BART);

print(L8TOA)
/*
// Mask the image collection
var L8EVImask = L8EVI.map(function(img) {
  var elev = elevGT500
  return img.updateMask(elev)
})

// Show the masked image colllection in the console
print(L8EVImask)

// Create EVI chart for flight box 
var series = ui.Chart.image.doySeriesByYear(
    L8EVImask, 'EVI', BART, ee.Reducer.mean(), 30);
    
// Display the chart
print(series);

// explicitly calculate the mean, store as a feature collection with no geometry, and export

// get the mean value for the region from each image
var ts = L8EVImask.map(function(image){
   var date = image.get('system:time_start');
   var mean = image.reduceRegion({
     reducer: ee.Reducer.mean(),
     geometry: BART,
     scale: 30
   });
   // and return a feature with 'null' geometry with properties (dictionary)  
   return ee.Feature(null, {'mean': mean.get('EVI'),
   //return ee.Feature(null, {'mean': mean,
                            'date': date})
});

// Show the time series feature collection in the console
print(ts);

// Export a .csv table of date, mean EVI for flight box
Export.table.toDrive({
  collection: ts,
  description: 'Mean_Landsat_EVI_2015-2018_XXXX',
  folder: 'EVI_Means',
  fileFormat: 'CSV'
});

// Display true color MODIS reflectance image  
//Map.addLayer(img.select(['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03']),
//         {gain: '0.1, 0.1, 0.1'}, 'MODIS bands 1/4/3');

*/
//Map.setCenter(-72.215200, 42.479365, 10); HARV
Map.setCenter(-71.271692, 44.036960, 10); BART
//Map.setCenter(-83.466381, 35.662826, 10); GRSM
//Map.setCenter(-84.470595, 31.242996, 10); JERC
//Map.setCenter(-66.853450, 17.986462, 10); GUAN
//Map.setCenter(-67.034225, 18.055970, 10); LAJA
//Map.setCenter(-89.508202, 46.207510, 10); UNDE
//Map.setCenter(-89.548808, 45.499229, 10); STEI_TREE
//Map.setCenter(-90.083689, 45.808634, 10); CHEQ
//Map.setCenter(-96.571292, 39.145820, 10); KONZ_KONA
//Map.setCenter(-95.192151, 39.022709, 10); UKFS
//Map.setCenter(-84.317737, 35.939079, 10); ORNL
//Map.setCenter(-80.531129, 37.381850, 10); MLBS
//Map.setCenter(-87.418400, 32.928635, 10); TALL
//Map.setCenter(-87.803750, 32.550871, 10); DELA
//Map.setCenter(-88.193000, 31.832444, 10); LENO
//Map.setCenter(-99.176052, 47.160925, 10); WOOD_DCFS
//Map.setCenter(-100.913440, 46.785939, 10); NOGP
//Map.setCenter(-97.608486, 33.369779, 10); CLBJ
//Map.setCenter(-99.100272, 35.363809, 10); OAES
//Map.setCenter(-110.520610, 44.919800, 10); YELL
//Map.setCenter(-105.567868, 40.027910, 10); NIWO
//Map.setCenter(-109.38827, 38.24833, 10); MOAB
//Map.setCenter(-106.84254, 32.59068, 10); JORN
//Map.setCenter(-110.880952, 31.840961, 10); SRER
//Map.setCenter(-112.491058, 40.180352, 10); ONAQ
//Map.setCenter(-119.731645, 37.090187, 10); SJER
//Map.setCenter(-119.006613, 37.01607, 10); TEAK
//Map.setCenter(-149.397424, 68.592491, 10); TOOL
//Map.setCenter(-156.620476, 71.230223, 10); BARR
//Map.setCenter(-147.508713, 65.184184, 10); BONA
//Map.setCenter(-145.768081, 63.858357, 10); DEJU
//Map.setCenter(-149.21334, 63.879976, 10); HEAL
//Map.setCenter(-104.734865, 40.826763, 10); CPER
//Map.setCenter(-103.042, 40.479985, 10); STER
//Map.setCenter(-105.493606, 40.234759, 10); RMNP
//Map.setCenter(-121.951912, 45.820488, 10); WREF
//Map.setCenter(-122.311358, 45.762317, 10); ABBY
//Map.setCenter(-155.271043, 19.557675, 10); OLAA

// Display an EVI 8-day composite image date
var dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI')
                  .filterDate('2015-08-09', '2015-08-17');
var colorized = dataset.select('EVI');
var colorizedVis = {
  min: 0.0,
  max: 1.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};
Map.addLayer(colorized, colorizedVis, 'Colorized');
print(colorized)
var vizParams = {
  bands: ['B5', 'B4', 'B3'],
  min: 0,
  max: 0.5,
  gamma: [0.95, 1.1, 1]};

Map.addLayer(L8TOA, vizParams)
//Map.addLayer(L8EVI,colorizedVis, 'L8EVI')
//Map.addLayer(L8EVImask, colorizedVis, 'L8EVImask')
//Display site table 
Map.addLayer(BART)
