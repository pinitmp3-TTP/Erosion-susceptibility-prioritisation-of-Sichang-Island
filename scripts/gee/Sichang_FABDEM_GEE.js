// =====================================================================
// Sichang Island: FABDEM (Copernicus GLO-30 with forests and buildings
// removed; Hawker et al., 2022) on the analysis grid, for a DEM sensitivity test
// Grid: EPSG:32647, 30.66 m, same origin as dem.tif (261 x 130 cells)
// Export: Google Drive, folder GEE_Sichang
// =====================================================================
var fab = ee.ImageCollection('projects/sat-io/open-datasets/FABDEM').mosaic().rename('fabdem');
var glo = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').mosaic().rename('glo30');

var transform = [30.66, 0, 694739.4368685539, 0, -30.66, 1457976.0128215954];
var region = ee.Geometry.Rectangle([694739.4368685539, 1449973.1843870878,
                                    698725.5199968526, 1457976.0128215954], 'EPSG:32647', false);

var both = fab.addBands(glo).toFloat();
Map.centerObject(region, 13);
Map.addLayer(fab.subtract(glo), {min: -15, max: 2, palette: ['b2182b', 'f7f7f7', '2166ac']},
             'FABDEM - GLO30 (m)');

Export.image.toDrive({
  image: both,
  description: 'Sichang_FABDEM_grid',
  folder: 'GEE_Sichang',
  fileNamePrefix: 'Sichang_FABDEM_grid',
  region: region,
  crs: 'EPSG:32647',
  crsTransform: transform,
  maxPixels: 1e9
});
