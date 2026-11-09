// =====================================================================
// Sichang Island: Sentinel-2 bare-soil frequency, 2019-2024
// Bare observation = NDVI < threshold AND bare soil index (BSI) > 0
// Frequency = bare observations / clear observations (Cloud Score+ cs_cdf >= 0.6)
// Computed at 10 m, clipped to the island at 10 m, averaged to the 30 m export grid (EPSG:32647).
// Output bands: freq15, freq20, freq25 (NDVI thresholds 0.15/0.20/0.25), n_valid
// Export: Google Drive, folder GEE_Sichang
// =====================================================================
var island = ee.FeatureCollection('projects/sichang-island/assets/Sichang_Island');
var region = ee.Geometry.Rectangle([692880, 1449960, 699480, 1460010], 'EPSG:32647', false);

var csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(island)
  .filterDate('2019-01-01', '2025-01-01')
  .linkCollection(csp, ['cs_cdf']);

var THR = [0.15, 0.20, 0.25];

var perImage = s2.map(function (img) {
  var clear = img.select('cs_cdf').gte(0.6);
  var r = img.select(['B2', 'B4', 'B8', 'B11']).divide(10000);
  var ndvi = r.normalizedDifference(['B8', 'B4']);
  var bsi = r.expression('((S + R) - (N + B)) / ((S + R) + (N + B))',
                         {S: r.select('B11'), R: r.select('B4'), N: r.select('B8'), B: r.select('B2')});
  var bands = THR.map(function (t) {
    return ndvi.lt(t).and(bsi.gt(0)).rename('b' + Math.round(t * 100));
  });
  return ee.Image.cat(bands).addBands(ee.Image(1).rename('valid'))
           .updateMask(clear).unmask(0, false)
           .multiply(clear)        // count only clear observations
           .addBands(clear.rename('nclear'));
});

var sums = perImage.sum();
var n = sums.select('nclear');
var freq = sums.select(['b15', 'b20', 'b25']).divide(n.max(1))
  .rename(['freq15', 'freq20', 'freq25'])
  .addBands(n.rename('n_valid'))
  .updateMask(n.gt(0));

// Clip to the island at 10 m (so sea and beach pixels are excluded), then average to 30 m
var island10 = ee.Image(1).clip(island).reproject({crs: 'EPSG:32647', scale: 10});
var freq10 = freq.updateMask(island10.mask())
  .reproject({crs: 'EPSG:32647', crsTransform: [10, 0, 692880, 0, -10, 1460010]});
var freq30 = freq10.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 64})
  .reproject({crs: 'EPSG:32647', crsTransform: [30, 0, 692880, 0, -30, 1460010]})
  .toFloat();

Map.centerObject(island, 14);
Map.addLayer(freq30.select('freq15'), {min: 0, max: 0.3, palette: ['ffffff', 'fdae6b', '7f2704']}, 'Bare frequency (NDVI < 0.15)');
print('Images used', s2.size());

Export.image.toDrive({
  image: freq30,
  description: 'Sichang_bare_frequency_30m',
  folder: 'GEE_Sichang',
  fileNamePrefix: 'Sichang_bare_frequency_30m',
  region: region,
  crs: 'EPSG:32647',
  crsTransform: [30, 0, 692880, 0, -30, 1460010],
  maxPixels: 1e9
});
