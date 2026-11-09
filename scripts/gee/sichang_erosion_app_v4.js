// ================================================================
// Ko Si Chang: erosion-susceptibility priority (v4)
// Two-block framework: RUSLE supply + Index of Connectivity (IC)
// 41 units = 29 complete sea-outlet basins (B) + 12 coastal faces (F)
//
// Assets required (upload first to projects/sichang-island/assets/):
//   sichang_units_v4    <- sichang_units_v4.zip  (shapefile, 41 units)
//   sichang_rasters_v4  <- sichang_rasters_v4.tif (b1 = soil loss A, b2 = IC, b3 = bare-soil frequency)
//   Sichang_Island      (already uploaded for v3)
// ================================================================
var ROOT = 'projects/sichang-island/assets/';
var UNITS = ee.FeatureCollection(ROOT + 'sichang_units_v4');
var RAS = ee.Image(ROOT + 'sichang_rasters_v4').select(['b1', 'b2', 'b3'], ['A', 'IC', 'bare']);
var ISLAND = ee.FeatureCollection(ROOT + 'Sichang_Island');
var RAS_M = RAS.updateMask(RAS.gt(-9999));

var PAL = ['1a9641', 'a6d96a', 'ffffbf', 'fdae61', 'd7191c'];
var PRIO_NAMES = ['Very low', 'Low', 'Moderate', 'High', 'Very high (top quintile)'];
var LOSS_NAMES = ['< 12.5', '12.5–31.25', '31.25–93.75', '93.75–125', '> 125'];
var IC_PAL = ['440154', '3b528b', '21918c', '5ec962', 'fde725'];
var BARE_PAL = ['ffffff', 'fee6ce', 'fdae6b', 'e6550d', '7f2704'];

function lossClass(img) {
  return img.gt(12.5).add(img.gt(31.25)).add(img.gt(93.75)).add(img.gt(125)).add(1)
            .updateMask(img.mask()).rename('class');
}

// ---------------- Map layers ----------------
Map.setOptions('SATELLITE');
Map.style().set('cursor', 'crosshair');

var prioImg = ee.Image().byte().paint(UNITS, 'priority');
var outline = ee.Image().byte().paint({featureCollection: UNITS, color: 1, width: 1.5});

var layers = {
  prio:  ui.Map.Layer(prioImg, {min: 1, max: 5, palette: PAL}, 'Priority class', true, 0.65),
  loss:  ui.Map.Layer(lossClass(RAS_M.select('A')), {min: 1, max: 5, palette: PAL}, 'Supply: RUSLE soil loss', false, 0.85),
  ic:    ui.Map.Layer(RAS_M.select('IC'), {min: -8, max: -2, palette: IC_PAL}, 'Transfer: Index of Connectivity', false, 0.85),
  bare:  ui.Map.Layer(RAS_M.select('bare'), {min: 0, max: 0.3, palette: BARE_PAL}, 'Bare-soil frequency 2019–2024', false, 0.85),
  line:  ui.Map.Layer(outline, {palette: ['000000']}, 'Unit boundaries', true),
  coast: ui.Map.Layer(ee.Image().byte().paint({featureCollection: ISLAND, color: 1, width: 2}),
                      {palette: ['FFFFFF']}, 'Island coastline', true)
};
['loss', 'ic', 'bare', 'prio', 'line', 'coast'].forEach(function (k) { Map.layers().add(layers[k]); });
var hiLayer = null;

// ---------------- Helpers ----------------
function swatchLegend(title, names, pal) {
  var p = ui.Panel({style: {margin: '6px 0'}});
  p.add(ui.Label(title, {fontWeight: 'bold', fontSize: '13px', margin: '4px 0'}));
  names.forEach(function (n, i) {
    p.add(ui.Panel([
      ui.Label('', {backgroundColor: '#' + pal[i], padding: '8px', margin: '2px 6px 2px 0'}),
      ui.Label(n, {fontSize: '12px', margin: '4px 0'})
    ], ui.Panel.Layout.Flow('horizontal')));
  });
  return p;
}
function rampLegend(title, pal, lo, hi) {
  var bar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '200x12', format: 'png', min: 0, max: 1, palette: pal},
    style: {stretch: 'horizontal', margin: '2px 8px', maxHeight: '14px'}
  });
  return ui.Panel([
    ui.Label(title, {fontWeight: 'bold', fontSize: '13px', margin: '6px 0 2px 0'}), bar,
    ui.Panel([ui.Label(lo, {fontSize: '11px', margin: '0 8px'}),
              ui.Label(hi, {fontSize: '11px', margin: '0 8px', textAlign: 'right', stretch: 'horizontal'})],
             ui.Panel.Layout.Flow('horizontal'))
  ]);
}
function row(k, v) {
  return ui.Panel([
    ui.Label(k, {width: '185px', fontSize: '12px', margin: '1px 0'}),
    ui.Label(String(v), {fontSize: '12px', fontWeight: 'bold', margin: '1px 0'})
  ], ui.Panel.Layout.Flow('horizontal'));
}
function fmt(x, d) { return (x === null || x === undefined || isNaN(x)) ? '–' : Number(x).toFixed(d); }
function toggle(label, layer) {
  return ui.Checkbox({label: label, value: layer.getShown(), onChange: function (c) { layer.setShown(c); }});
}
function head(t) { return ui.Label(t, {fontWeight: 'bold', fontSize: '13px', margin: '10px 0 2px 0'}); }

// ---------------- Side panel ----------------
var panel = ui.Panel({style: {width: '380px', padding: '8px'}});
panel.add(ui.Label('Ko Si Chang: erosion-susceptibility priority', {fontSize: '18px', fontWeight: 'bold'}));
panel.add(ui.Label('Two-block ranking of 41 units (29 complete sea-outlet basins, B; 12 coastal faces, F). ' +
                   'Supply = RUSLE soil loss per erodible area; transfer = Index of Connectivity to the sea. ' +
                   'Click the map to inspect a unit.', {fontSize: '12px'}));

panel.add(head('Layers'));
panel.add(toggle('Priority class', layers.prio));
panel.add(toggle('Supply: RUSLE soil loss (t ha⁻¹ yr⁻¹)', layers.loss));
panel.add(toggle('Transfer: Index of Connectivity', layers.ic));
panel.add(toggle('Evidence: bare-soil frequency 2019–2024', layers.bare));
panel.add(toggle('Unit boundaries', layers.line));
panel.add(toggle('Island coastline', layers.coast));

panel.add(swatchLegend('Priority class (compound score quintiles)', PRIO_NAMES, PAL));
panel.add(swatchLegend('Soil loss (t ha⁻¹ yr⁻¹; LDD classes)', LOSS_NAMES, PAL));
panel.add(rampLegend('Index of Connectivity', IC_PAL, '−8 (poorly connected)', '−2 (well connected)'));
panel.add(rampLegend('Bare-soil frequency (built-up areas left blank)', BARE_PAL, '0', '≥ 0.3'));

var info = ui.Panel({style: {margin: '8px 0', padding: '6px', border: '1px solid #cccccc'}});
info.add(ui.Label('Click on the map to inspect a unit.', {color: '666666', fontSize: '12px'}));
panel.add(head('Selected unit'));
panel.add(info);

var topUnits = UNITS.sort('CS', false).limit(12);
panel.add(ui.Chart.feature.byFeature({features: topUnits, xProperty: 'label', yProperties: ['CS']})
  .setChartType('BarChart')
  .setOptions({title: 'Top 12 units by compound score', legend: {position: 'none'},
               hAxis: {title: 'Compound score (z)', minValue: 0},
               vAxis: {textStyle: {fontSize: 11}},
               chartArea: {left: 45, top: 30, width: '75%', height: '80%'},
               colors: ['#d7191c'], height: 330}));

panel.add(ui.Label(
  'How to read this map: the ranking shows relative susceptibility from topography and cover. ' +
  'The top-priority units are mostly intact, steep forest, so they mark where clearing would matter most. ' +
  'Present-day bare soil and post-rain turbidity were higher outside these units. ' +
  'Indicative only; not field-validated.',
  {fontSize: '11px', color: '666666', margin: '10px 0 0 0'}));

ui.root.insert(0, panel);
Map.setCenter(100.8094, 13.1527, 14);

// ---------------- Click inspector ----------------
Map.onClick(function (c) {
  var pt = ee.Geometry.Point([c.lon, c.lat]);
  info.clear();
  info.add(ui.Label('Loading…', {color: '666666', fontSize: '12px'}));
  var hit = UNITS.filterBounds(pt);
  var props = ee.Algorithms.If(hit.size().gt(0), ee.Feature(hit.first()).toDictionary(), ee.Dictionary({}));
  var pix = RAS_M.reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 30});

  if (hiLayer) { Map.layers().remove(hiLayer); }
  hiLayer = ui.Map.Layer(hit.style({color: '00FFFF', fillColor: '00000000', width: 3}), {}, 'Selected unit');
  Map.layers().add(hiLayer);

  ee.Dictionary({props: props, pix: pix}).evaluate(function (r, err) {
    info.clear();
    if (err) { info.add(ui.Label('Error: ' + err, {color: 'red', fontSize: '12px'})); return; }
    var p = r.props, x = r.pix;
    if (!p || Object.keys(p).length === 0) {
      info.add(ui.Label('Outside the ranked units (sea).', {fontSize: '12px'}));
    } else {
      var kind = p.kind === 'basin' ? 'complete sea-outlet basin' : 'coastal face';
      info.add(ui.Label('Unit ' + p.label + ' (' + kind + ')', {fontWeight: 'bold', fontSize: '14px'}));
      info.add(row('Priority class', p.prio_lab));
      info.add(row('Rank (of 41)', fmt(p.rank, 0)));
      info.add(row('P(top quintile), Monte Carlo', fmt(p.P_top, 2)));
      info.add(row('Dominant driver', p.driver));
      info.add(row('Area (ha)', fmt(p.area_ha, 1)));
      info.add(row('Mean slope (°)', fmt(p.slope_deg, 1)));
      info.add(row('Supply (t ha⁻¹ yr⁻¹)', fmt(p.supply, 1)));
      info.add(row('Mean IC', fmt(p.IC, 2)));
      info.add(row('z(supply) / z(transfer)', fmt(p.z_supply, 2) + ' / ' + fmt(p.z_trans, 2)));
      info.add(row('Compound score', fmt(p.CS, 2)));
      info.add(ui.Label('Independent evidence', {fontSize: '12px', fontWeight: 'bold', margin: '6px 0 1px 0'}));
      info.add(row('Bare-soil frequency', fmt(p.bare_freq, 3)));
      info.add(row('Turbidity response after rain (Δln)', fmt(p.turb_resp, 2)));
      if (p.kind === 'basin') {
        info.add(ui.Label('Basin morphometry', {fontSize: '12px', fontWeight: 'bold', margin: '6px 0 1px 0'}));
        info.add(row('Drainage density Dd (km km⁻²)', fmt(p.Dd, 2)));
        info.add(row('Relief H (m)', fmt(p.H, 0)));
        info.add(row('Relief ratio Rh', fmt(p.Rh, 3)));
        info.add(row('Hypsometric integral HI', fmt(p.HI, 2)));
        info.add(row('Form factor Ff / circularity Rc', fmt(p.Ff, 2) + ' / ' + fmt(p.Rc, 2)));
      }
    }
    if (x && x.A !== null && x.A !== undefined) {
      info.add(ui.Label('At clicked point: soil loss ' + fmt(x.A, 1) + ' t ha⁻¹ yr⁻¹, IC ' + fmt(x.IC, 2) +
                        ', bare-soil frequency ' + fmt(x.bare, 3),
                        {fontSize: '11px', color: '444444', margin: '6px 0 0 0'}));
    }
  });
});
