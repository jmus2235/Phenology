// Earth Engine script to produce time series of 8-day composite mean EVI estimates from MODIS and EVI2 from VIIRS within 
//  ROIs defined by NEON AOP flight boxes, TOS boundaries, vegetation types, etc. 
//  Region of Interests (ROI) asset tables must be uploaded from shapefiles into GEE and imported into script. 
//  Example below shows variable siteID defined as the JORN flight box. 
//  Output from this script are 3 three-column tables: a 16-day composite from MODIS (Terra) and a 16-day composite MODIS (Aqua)
//  to be combined into a single table 8-day composite table; and a third table from VIIRS (already 8-day composite). 
//  These tables are used as input into R scripts "MODIS_EVI_processing_v1-2.R" and "VIIRS_EVI_processing_v1-2.R", respectively. 

var siteID = JORN
var EVIstartDate = ('2003-01-01')
var EVIendDate = ('2019-12-31')
var VIIRStartDate = ('2012-01-17')

/* Create a function that select images on basis of QA information (stored in QA band).
Return only images according to specified bits. 

  Function Arguments:
*   image - The QA Image to get bits from.
*   start - The first bit position, 0-based.
*   end   - The last bit position, inclusive.
*   name  - A name for the output image.
*/

var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
      pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
};

// Define a function to mask out poor quality pixels (bits 0 and 1) in MODIS
var mask= function(image) {
//Select the QA band.
  var QA = image.select('DetailedQA');
  var cloud = getQABits(QA, 0, 1, 'VI quality')
                        //.expression("b(0) == 2 || b(0) == 3");
                        .expression("b(0) == 1 || b(0) == 2 || b(0) == 3");
 // Return an image masking out poor quality (cloudy) areas
return image.updateMask(cloud.not())//.internalCloud.eq(0));
};

// Define a function to mask out poor quality pixels (bits 0 and 1) in VIIRS
var maskVIIRS= function(image) {
//Select the QA band.
  var QA = image.select('VI_Quality');
  var cloudVIIRS = getQABits(QA, 0, 1, 'MODLAND_QA')
                        //.expression("b(0) == 2 || b(0) == 3");
                        .expression("b(0) == 1 || b(0) == 2 || b(0) == 3");
 // Return an image masking out poor quality (cloudy) areas
return image.updateMask(cloudVIIRS.not())//.internalCloud.eq(0));
};

//EVI
//load the MODIS TERRA image collection
var modisEVI = ee.ImageCollection('MODIS/006/MOD13Q1') //MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m
    .filter(ee.Filter.dayOfYear(0, 365))
    .filterDate(EVIstartDate, EVIendDate)
    .filterBounds(siteID)
    .map(mask)

var series = ui.Chart.image.doySeriesByYear(
    modisEVI, 'EVI', siteID, ee.Reducer.mean(), 250);

print(modisEVI)    
// Display the chart
print(series);

// explicitly calculate the mean, store as a feature collection with no geometry, and export

// get the mean value for the region from each image
var ts = modisEVI.map(function(image){
   var date = image.get('system:time_start');
   var mean = image.reduceRegion({
     reducer: ee.Reducer.mean(),
     geometry: siteID,
     scale: 250
   });
   // and return a feature with 'null' geometry with properties (dictionary)  
   return ee.Feature(null, {'mean': mean.get('EVI'),
   //return ee.Feature(null, {'mean': mean,
                            'date': date})
});

print(ts);

// Export a .csv table of date, mean EVI for flight box
Export.table.toDrive({
  collection: ts,
  description: 'Mean_TERRA_EVI_2003-2019_MOD13Q1_0-bit_XXXX',
  folder: 'EVI_Means',
  fileFormat: 'CSV'
});

//load the MODIS AQUA image collection
var modisEVI = ee.ImageCollection('MODIS/006/MYD13Q1') //MYD13Q1.006 AQUA Vegetation Indices 16-Day Global 250m
    .filter(ee.Filter.dayOfYear(0, 365))
    .filterDate(EVIstartDate, EVIendDate)               
    .filterBounds(siteID)                       
    .map(mask)
    
var series = ui.Chart.image.doySeriesByYear(
    modisEVI, 'EVI', siteID, ee.Reducer.mean(), 250);
    
// Display the chart
print(series);

// explicitly calculate the mean, store as a feature collection with no geometry, and export

// get the mean value for the region from each image
var ts = modisEVI.map(function(image){
   var date = image.get('system:time_start');
   var mean = image.reduceRegion({
     reducer: ee.Reducer.mean(),
     geometry: siteID,
     scale: 250
   });
   // and return a feature with 'null' geometry with properties (dictionary)  
   return ee.Feature(null, {'mean': mean.get('EVI'),
                            'date': date})
});

print(ts);

// Export a .csv table of date, mean EVI for flight box
Export.table.toDrive({
  collection: ts,
  description: 'Mean_AQUA_EVI_2003-2019_MYD13Q1_0-bit_XXXX',
  folder: 'EVI_Means',
  fileFormat: 'CSV'
});

//load VIIRS image collection
var viirsEVI = ee.ImageCollection('NOAA/VIIRS/001/VNP13A1') //VIIRS Vegetation Indices 16-Day Global 500m
    .filter(ee.Filter.dayOfYear(0, 365))
    .filterDate(VIIRStartDate, EVIendDate)               
    .filterBounds(siteID)                       
    .map(maskVIIRS)
    
var series = ui.Chart.image.doySeriesByYear(
    viirsEVI, 'EVI2', siteID, ee.Reducer.mean(), 500);
    
// Display the chart
print(series);

// explicitly calculate the mean, store as a feature collection with no geometry, and export

// get the mean value for the region from each image
var ts2 = viirsEVI.map(function(image){
   var date = image.get('system:time_start');
   var mean = image.reduceRegion({
     reducer: ee.Reducer.mean(),
     geometry: siteID,
     scale: 500
   });
   // and return a feature with 'null' geometry with properties (dictionary)  
   return ee.Feature(null, {'mean': mean.get('EVI2'),
                            'date': date})
});

print(ts2);

// Export a .csv table of date, mean EVI for flight box
Export.table.toDrive({
  collection: ts2,
  description: 'Mean_VIIRS_EVI2_2012-2019_VNP13A1_0-bit_XXXX',
  folder: 'EVI_Means',
  fileFormat: 'CSV'
});


//Map.addLayer(BART_DB)
//Map.setCenter(-71.271877, 44.036956, 10); BART 
//Map.addLayer(DSNY_TOS)
//Map.setCenter(-81.431683, 28.088751, 10); DSNY 
