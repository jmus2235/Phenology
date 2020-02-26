// Script that calculates number of cloudy MODIS pixels for each day within an ROI. 
// This is a beta script that produces a table requiring significant manual re-processing in Excel.
// Final table appears in the following format, which is then meant to be input into a pivot table:
// YEAR	MONTH	DAY	GROUP 0	COUNT 0	GROUP 1	COUNT 1	TOTAL	CLOUD FREE FRACTION	GREATER THAN 90% CLOUD FREE
// 2002	1	    1	        0	   1227	      1	     21	 1248	        0.983173077	                          1
// 2002	1	    2	        0	   1248			               1248	                  1	                          1
// 2002	1	    3	        0	   1248			               1248	                  1                          	1
// ...

// Modis Cloud Masking 
//
// Calculate how frequently a location is labeled as clear (i.e. non-cloudy)
// according to the "internal cloud algorithm flag" of the MODIS "state 1km"
// QA band.

/*
 * A function that returns an image containing just the specified QA bits.
 *
 * Args:
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

// A function to mask out pixels that did not have observations.
var maskEmptyPixels = function(image) {

// Find pixels that had observations.
  var withObs = image.select('num_observations_1km').gt(0);
    //Selects bands from an image.
    //Returns an image with the selected bands.
  
  return image.updateMask(withObs); 
    //Updates an image's mask at all positions where the existing mask is not zero. 
    //The output image retains the metadata and footprint of the input image.
};

// A function to mask out cloudy pixels.
var maskClouds = function(image) {
  // Select the QA band.
  var QA = image.select('state_1km');
  // Get the internal_cloud_algorithm_flag bit.
  var cloud = getQABits(QA, 10, 10, 'internal_cloud_algorithm_flag').neq(0).rename('cloud')
  return image.addBands(cloud)
};

// Start with a MODIS TERRA image collection for an X day period
//var collection = ee.ImageCollection('MODIS/006/MOD09GA')
//                   .filterDate('2001-01-01', '2009-12-31');

// Start with an MODIS AQUA image collection for an X day period
var collection = ee.ImageCollection('MODIS/006/MOD09GA')
//                   .filterDate('2002-01-01', '2009-12-31');
                   .filterDate('2010-01-01', '2017-12-31');

// Initialize an empty Dictionary
//var fractionInit = ee.Dictionary();
                   
var cloudy = collection.map(maskClouds)
var result = cloudy.map(function(image) {
  var cloud = image.select('cloud')
  return ee.Feature(null, cloud.addBands(cloud).reduceRegion({
    reducer: ee.Reducer.count().group(),
    geometry: SMAP_MB.geometry()
  }))
})
print(result)

//var result2 = collection.map(function(f) { 
//   return ee.Feature(f.geometry(), f.toDictionary().combine({count: 0, group: 0}, false)
//})

//var cloudfractions=ee.Dictionary(result);
//print('Dictionary of fractions', result);

//Export.table.toDrive(ee.FeatureCollection(result), 'SiteCloudCounts');

// Make a feature without geometry and set the properties to the dictionary of results
//var feature = ee.Feature(null, result);

// Wrap the Feature in a FeatureCollection for export.
//var featureCollection = ee.FeatureCollection([result]).flatten();

//Export the feature collection
Export.table.toDrive({
  collection: result,
  description: 'SiteCloudCounts',
  folder: 'Cloud_Fractions',
  fileNamePrefix: 'SiteCloudCounts',
  fileFormat: 'CSV'
});


// Call the MODIS reflectance data for this date (JM addition)
var img = ee.Image('MOD09GA/MOD09GA_005_2016_09_14');

//Map.addLayer(
//    clearObsCount.toFloat().divide(totalObsCount),
//    {min: 0, max: 1},
//    'ratio of clear to total observations'
//  );

// Display true color MODIS reflectance image  
Map.addLayer(img.select(['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03']),
         {gain: '0.1, 0.1, 0.1'}, 'MODIS bands 1/4/3');

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
//Map.setCenter(-122.023675, 37.279722, 10); SCMS
//Map.setCenter(-109.583784, 38.089812, 10); CTW1
//Map.setCenter(-112.064353, 34.249008, 10); CTW2
//Map.setCenter(-112.847911, 34.919718, 10); CTW3
//Map.setCenter(-114.469337, 32.858003, 10); CTW4
//Map.setCenter(-90.869677, 29.322723, 10); DLTX
Map.setCenter(-72.123258, 42.519709, 10); SMAP_MA

//Display site table 
Map.addLayer(SMAP_MA)
