// =====================================================================
// Sichang Island: nearshore turbidity at sub-basin outlets after rain
// Sentinel-2 L2A (harmonised) 2019-2025, Cloud Score+, GPM IMERG antecedent rain
// Output: one row per image x point (outlet and paired offshore control)
// Run, then click "Run" on the task in the Tasks tab to export the CSV.
// =====================================================================

// ---- 1. Points: one outlet per unit (60-150 m seaward of the coast)
//         and a paired control 1 km further offshore along the same normal
var pts = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([100.805625, 13.142025]), {ucode: 'U01', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.796734, 13.144426]), {ucode: 'U01', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.80264, 13.179925]), {ucode: 'U02', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.802706, 13.188964]), {ucode: 'U02', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.801304, 13.172985]), {ucode: 'U03', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.793284, 13.168522]), {ucode: 'U03', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.809511, 13.170433]), {ucode: 'U04', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.817466, 13.165858]), {ucode: 'U04', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.819115, 13.151264]), {ucode: 'U05', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.819182, 13.160302]), {ucode: 'U05', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.814889, 13.155836]), {ucode: 'U06', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.824111, 13.155771]), {ucode: 'U06', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.816153, 13.153483]), {ucode: 'U07', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.820822, 13.161278]), {ucode: 'U07', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.820332, 13.139446]), {ucode: 'U08', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.829554, 13.139381]), {ucode: 'U08', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.803405, 13.150781]), {ucode: 'U09', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.79545, 13.155356]), {ucode: 'U09', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.816666, 13.152583]), {ucode: 'U10', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.823235, 13.158927]), {ucode: 'U10', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.816316, 13.131715]), {ucode: 'U11', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.825537, 13.13165]), {ucode: 'U11', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.819069, 13.147492]), {ucode: 'U12', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.828291, 13.147427]), {ucode: 'U12', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.801015, 13.172138]), {ucode: 'U13', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.793061, 13.176713]), {ucode: 'U13', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.811755, 13.126758]), {ucode: 'U14', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.820977, 13.126693]), {ucode: 'U14', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.80362, 13.155084]), {ucode: 'U15', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.794397, 13.155149]), {ucode: 'U15', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.80616, 13.143426]), {ucode: 'U16', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.796937, 13.143491]), {ucode: 'U16', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.815874, 13.154039]), {ucode: 'U17', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.820543, 13.161834]), {ucode: 'U17', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.813657, 13.15736]), {ucode: 'U18', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.821678, 13.161823]), {ucode: 'U18', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.806419, 13.134683]), {ucode: 'U19', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.80175, 13.126888]), {ucode: 'U19', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.806494, 13.144937]), {ucode: 'U20', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.801825, 13.137142]), {ucode: 'U20', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.812242, 13.127779]), {ucode: 'U21', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.818716, 13.121342]), {ucode: 'U21', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.817457, 13.133092]), {ucode: 'U22', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.826679, 13.133027]), {ucode: 'U22', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.806194, 13.145514]), {ucode: 'U23', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.798174, 13.141051]), {ucode: 'U23', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.805871, 13.142597]), {ucode: 'U24', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.796649, 13.142662]), {ucode: 'U24', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.80874, 13.165857]), {ucode: 'U25', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.817964, 13.165791]), {ucode: 'U25', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.811422, 13.161275]), {ucode: 'U26', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.817991, 13.16762]), {ucode: 'U26', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.806466, 13.144108]), {ucode: 'U27', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.798513, 13.148683]), {ucode: 'U27', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.819049, 13.144721]), {ucode: 'U28', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.828271, 13.144656]), {ucode: 'U28', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.809228, 13.167793]), {ucode: 'U29', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.818451, 13.167728]), {ucode: 'U29', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.809448, 13.164616]), {ucode: 'U30', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.814117, 13.172411]), {ucode: 'U30', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.815313, 13.154598]), {ucode: 'U31', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.819982, 13.162392]), {ucode: 'U31', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.819247, 13.145828]), {ucode: 'U32', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.82847, 13.145763]), {ucode: 'U32', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.81333, 13.160366]), {ucode: 'U33', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.819899, 13.166711]), {ucode: 'U33', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.820499, 13.149699]), {ucode: 'U34', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.829722, 13.149634]), {ucode: 'U34', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.810558, 13.161837]), {ucode: 'U35', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.815228, 13.169632]), {ucode: 'U35', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.813677, 13.157508]), {ucode: 'U36', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.8229, 13.157442]), {ucode: 'U36', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.8098, 13.174036]), {ucode: 'U37', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.814469, 13.181831]), {ucode: 'U37', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.810295, 13.161837]), {ucode: 'U38', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.816864, 13.168182]), {ucode: 'U38', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.816277, 13.131842]), {ucode: 'U39', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.82083, 13.123982]), {ucode: 'U39', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.804331, 13.179173]), {ucode: 'U40', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.804398, 13.188211]), {ucode: 'U40', kind: 'control'}),
  ee.Feature(ee.Geometry.Point([100.808935, 13.16641]), {ucode: 'U41', kind: 'outlet'}),
  ee.Feature(ee.Geometry.Point([100.818158, 13.166344]), {ucode: 'U41', kind: 'control'})
]);

var BUF      = 120;              // sampling radius around each point (m)
var START    = '2019-01-01';
var END      = '2025-11-01';
var CS_MIN   = 0.60;             // Cloud Score+ cs_cdf threshold
var aoi      = pts.geometry().bounds().buffer(1500);

// ---- 2. Sentinel-2 with Cloud Score+
var csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi).filterDate(START, END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .linkCollection(csp, ['cs_cdf']);

// ---- 3. Water-pixel masks and turbidity
function prep(img) {
  var r = img.select(['B2','B3','B4','B8','B11']).divide(10000);
  var clear = img.select('cs_cdf').gte(CS_MIN);
  var mndwi = r.normalizedDifference(['B3','B11']);
  var water = mndwi.gt(0.2)
      .and(r.select('B8').lt(0.06))      // drops ships, wakes, land
      .and(r.select('B11').lt(0.02));    // drops sun glint
  // remove mixed pixels at the land-water edge (2 pixels)
  water = water.focalMin({radius: 20, units: 'meters'}).selfMask();
  var rho = r.select('B4');
  // Nechad et al. (2009) single-band turbidity, 665 nm (FNU)
  var T = rho.multiply(228.1).divide(ee.Image(1).subtract(rho.divide(0.1641))).rename('T');
  var ndti = r.normalizedDifference(['B4','B3']).rename('NDTI');
  return ee.Image(T.addBands(ndti).addBands(r.select('B4').rename('rho665'))
          .updateMask(water).updateMask(clear)
          .copyProperties(img, ['system:time_start','MGRS_TILE','SPACECRAFT_NAME']));
}

// ---- 4. Antecedent rain from GPM IMERG V07 half-hourly (mm), island area
//         (IMERG covers the sea; CHIRPS is land-only and may miss the island)
var imerg = ee.ImageCollection('NASA/GPM_L3/IMERG_V07').select('precipitation');
var island = ee.Geometry.Point([100.808, 13.155]).buffer(10000);
function rainSum(d, nDays) {      // rain in the nDays*24 h before the overpass
  return imerg.filterDate(d.advance(-nDays, 'day'), d).sum().multiply(0.5)
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: island, scale: 11132})
    .get('precipitation');
}

// ---- 5. Sample every image at every point
var rows = s2.map(prep).map(function(img) {
  var d = ee.Date(img.get('system:time_start'));
  var P1 = rainSum(d, 1), P3 = rainSum(d, 3), P7 = rainSum(d, 7);
  var red = ee.Reducer.median().combine(ee.Reducer.count(), '', true);
  return img.reduceRegions({collection: pts.map(function(f) { return f.buffer(BUF); }),
                             reducer: red, scale: 10})
    .map(function(f) {
      return f.setGeometry(null).set({date: d.format('YYYY-MM-dd'),
        tile: img.get('MGRS_TILE'), sat: img.get('SPACECRAFT_NAME'),
        P1: P1, P3: P3, P7: P7});
    });
}).flatten().filter(ee.Filter.gt('T_count', 5));

print('images in collection', s2.size());
print('example rows', rows.limit(5));

// ---- 6. Map check
Map.centerObject(aoi, 13);
Map.addLayer(s2.median().select(['B4','B3','B2']), {min: 0, max: 2500}, 'S2 median RGB', false);
Map.addLayer(pts.filter(ee.Filter.eq('kind','outlet')).map(function(f){return f.buffer(BUF);}), {color: 'red'}, 'outlet buffers');
Map.addLayer(pts.filter(ee.Filter.eq('kind','control')).map(function(f){return f.buffer(BUF);}), {color: 'blue'}, 'control buffers');

// ---- 7. Export
Export.table.toDrive({
  collection: rows,
  description: 'Sichang_turbidity_outlets',
  folder: 'GEE_Sichang',
  fileFormat: 'CSV',
  selectors: ['date','tile','sat','ucode','kind','T_median','T_count','NDTI_median','rho665_median','P1','P3','P7']
});
