# Erosion-susceptibility prioritisation of Sichang Island (Gulf of Thailand): data, scripts and remote evaluation

Data and code for:

> Tanachaichoksirikun, P.,  Sirikaew, U.,  Seeboonruang, U., Danpipattrakul, W.,  Chen, W. Hillslope supply or sediment connectivity? A two-block framework and open Earth Engine tool for prioritising erosion-susceptible units on a data-scarce tropical island. *CATENA* (submitted). [Update with the article DOI on publication.]

Interactive application: https://sichang-island.projects.earthengine.app/view/sichang-erosion

Contact: Pinit Tanachaichoksirikun, King Mongkut's Institute of Technology Ladkrabang (pinit.ta@kmitl.ac.th)

## Contents

```
inputs/           input rasters and tables on the analysis grid
scripts/python/   processing and analysis (run in the order below)
scripts/gee/      Google Earth Engine scripts and the application assets
outputs/          all outputs produced by the scripts
interpretation/   blind visual interpretation (random points and erosion-feature inventories)
```

### inputs/

| File | Description |
|---|---|
| `grid_land_mask.tif` | Analysis grid (UTM 47N, 30.66 m) and land mask (1 = land; 6,923 cells, 650.9 ha) |
| `dem.tif`, `slope.tif` | Copernicus GLO-30 DEM (conditioned) and slope (degrees) |
| `flowdir.tif`, `flowacc.tif` | D8 flow direction (ESRI coding) and flow accumulation (cells) |
| `Sichang_WorldCover2021_10m.tif` | ESA WorldCover 2021 v200 classes (10 m) |
| `Sichang_C_WorldCover_30m.tif` | 30 m grid used to aggregate WorldCover class fractions |
| `Sichang_C_NDVI_2021_30m.tif` | NDVI-based C factor (2021), used for a sensitivity variant |
| `S2_bare_frequency_2019_2024_30m.tif` | Sentinel-2 bare-soil frequency, 2019–2024 (export of `scripts/gee/Sichang_bare_frequency_GEE.js`): share of clear observations with NDVI below a threshold and BSI > 0. Bands `freq15`, `freq20`, `freq25` (NDVI thresholds 0.15, 0.20, 0.25) and `n_valid` (clear observations). Band 2 (`freq20`) is used in the analysis (`BARE_BAND` in `common.py`) |
| `soil_samples.csv` | Nine soil samples: coordinates, texture, bulk density, permeability and K for 0.5, 1.0 and 2.0% organic matter (Horpeancharoen et al., 2022) |
| `Sichang_turbidity_outlets.csv` | Sentinel-2 nearshore turbidity at unit outlets and offshore controls, 2019–2025 (export of `scripts/gee/Sichang_turbidity_GEE.js`) |

Columns of `Sichang_turbidity_outlets.csv`: `date`; `tile`; `sat`; `ucode` (unit code); `kind` (`outlet` or `control`); `T_median` (Nechad et al., 2010 turbidity, FNU, median within 120 m); `T_count` (valid water pixels); `NDTI_median`; `rho665_median`; `P1`, `P3`, `P7` (GPM IMERG V07 rainfall over the previous 1, 3 and 7 days, mm). Unit codes are linked to unit labels in `outputs/units_5ha_scored.gpkg` and `interpretation/random_points_key.csv`.

### scripts/python/ (Python 3; numpy, pandas, geopandas, rasterio, scipy, scikit-learn, matplotlib)

Run from `scripts/python/`:

```
python s1_units.py                        # trace every land cell to its sea outlet
python s1b_units.py 5                     # ranking units: basins >= 5 ha + coastal faces (also 2 and 10)
python s1b_units.py 2
python s1b_units.py 10
python s3_rebuild.py                      # C (LDD values from WorldCover fractions), K, LS (MB, DG, DG capped)
python s2_icvariants.py                   # Index of Connectivity, four weightings
python s4_table.py 5                      # unit table: supply, IC, slope, bare-soil frequency, morphometry
python s4_table.py 2
python s4_table.py 10
python s5_montecarlo.py                   # Monte Carlo (1000 draws) -> outputs/units_5ha_mc_rerun.csv
python s8_points.py                       # regenerate the random interpretation points (seeded)
python v4_results.py                      # scores, classes, null model, bare soil, turbidity, random points, kappa
python analysis/roc_features.py           # AUC of feature inventories against maps
python s6_inventory.py                    # agreement between the two inventories (adds to v4_results.json)
python make_figs_v4.py                    # Figs 2-6 of the article
python s9_export.py                       # GeoTIFF and GeoPackage exports
```

Notes:
- The Monte Carlo results reported in the article are in `outputs/units_5ha_mc.csv`. `s5_montecarlo.py` reproduces them to within 0.01 in P(top) (random draws differ); `v4_results.py` reads the original file.
- `s8_points.py` regenerates the 405 random points exactly; the regenerated key is identical to `interpretation/random_points_key.csv`.

### scripts/gee/

| File | Description |
|---|---|
| `Sichang_turbidity_GEE.js` | Sentinel-2 L2A + Cloud Score+ turbidity at outlets and controls, with GPM IMERG antecedent rainfall |
| `Sichang_bare_frequency_GEE.js` | Sentinel-2 L2A + Cloud Score+ bare-soil frequency, 2019–2024, computed at 10 m within the island outline and averaged to 30 m |
| `sichang_erosion_app_v4.js` | The public Earth Engine application |
| `app_assets/` | Assets uploaded for the application (`sichang_units_v4.zip`, `sichang_rasters_v4.tif`: b1 = soil loss, b2 = IC, b3 = bare-soil frequency) and the script that makes them |

### outputs/

| File | Description |
|---|---|
| `units_5ha_scored.gpkg` / `.csv` | The 41 ranking units with supply, IC (four variants), standardised scores, compound score `CS`, priority class `cls` (5 = top quintile), dominant driver, Monte Carlo P(top) and rank percentiles, bare-soil frequency, slope and morphometry (basins) |
| `units_5ha_table.csv`, `units_2ha_table.csv`, `units_10ha_table.csv` | Unit tables for the three minimum basin sizes; `units_*ha.gpkg/.shp` are the unit polygons |
| `units_5ha_mc.csv` | Monte Carlo results used in the article |
| `rasters/` | Soil loss `A_*` (t ha⁻¹ yr⁻¹; capped Desmet–Govers, uncapped Desmet–Govers, Moore–Burch LS), `K_OM1`, `C_LDD`, `LS_*`, `IC_*` (four weightings), built-up and water fractions |
| `turbidity_unit_response.csv` | Median Δln turbidity on dry and wet dates and the response (wet − dry) per unit |
| `randompts_AB.csv` | Parsed random-point interpretations with derived indicators (`D` disturbed, `bare10` ≥ 10% bare, `sealed`, `closed` closed forest, `E` erosion feature) |
| `roc_features_results.csv` | AUC with 95% bootstrap CI by inventory and map |
| `v4_results.json` | All statistics reported in the article |
| `figs/` | Figures 2–6 |
| `*.npz`, `*.npy` | Intermediate arrays on the analysis grid |

### interpretation/

The two interpreters had not visited the island, were not involved in model development and had no access to the ranking. They are anonymised as Interpreter 1 (A in the scripts) and Interpreter 2 (B). The original forms were in Thai; categorical values are kept as entered and translated below. Free-text notes are in Thai.

| File | Description |
|---|---|
| `random_points_key.csv` | 405 random points: `pid`, unit (`uid`, `label`, `ucode`), WGS 84 and UTM coordinates, and the unit's class, compound score and P(top) |
| `random_points_interpreter1.csv`, `random_points_interpreter2.csv` | Interpretations of the random points (10 m circle) |
| `erosion_features_interpreter1.gpkg/.csv`, `erosion_features_interpreter2.gpkg/.csv` | Erosion features mapped in the tile-by-tile sweep (UTM 47N) |
| `sweep_tiles_300m.kml` | The 108 sweep tiles (300 m × 300 m) |

Random-point fields: `pid`; `lon`, `lat`; `img` (image date); `cover` (dominant cover); `bare` (bare-soil fraction, roofs and roads excluded); `dist` (disturbance in or touching the circle); `eros` (erosion feature within 50 m); `where` (feature location); `old_img` (older image date); `seen` (feature on older image); `conf` (confidence); `notes`.
Feature fields: `feat_id`; `class`; `lon`, `lat`; `img`; `size` (m² or m); `seen`; `old_img`; `conf`; `notes` (for Interpreter 1, the sweep tile ID).

Codebook (Thai → English):

| Field | Thai | English |
|---|---|---|
| cover | 1 ป่าเรือนยอดปิด | 1 closed-canopy forest |
| | 2 ไม้ต้น/ไม้พุ่มโปร่ง | 2 open trees / shrubs |
| | 3 หญ้า | 3 grass |
| | 4 ดินเปลือยธรรมชาติ | 4 natural bare ground |
| | 5 ดินเปลือยจากคนทำ | 5 human-made bare ground |
| | 6 ผิวแข็ง (ถนน/หลังคา/ลาน) | 6 sealed (road / roof / yard) |
| | 8 น้ำ | 8 water (excluded) |
| bare | 0 / น้อยกว่า 10% / 10–50% / มากกว่า 50% | 0 / < 10% / 10–50% / > 50% |
| dist | ไม่มี | none |
| | ถนน/ทางเดิน | road / track |
| | บ่อหิน/บ่อขุด | quarry / pit |
| | ก่อสร้าง/ถางใหม่ | construction / new clearing |
| | รอยไฟไหม้ | fire scar |
| eros, class | ไม่มี | none |
| | ร่องริ้ว/รอยน้ำพัดผิวดิน | rills / sheet-wash marks |
| | ร่องลึก/ร่องน้ำกัดลึก | gully |
| | รอยดินถล่ม | slope-failure scar |
| | ลาดตัดที่กำลังพัง | eroding road cut |
| where | ไม่มีร่องรอย / ในวงกลม / ภายใน 50 ม. | no feature / in the circle / within 50 m |
| seen | ไม่มีร่องรอย / มีแล้วในภาพเก่า / ภาพเก่ายังไม่มี | no feature / present on older image / not on older image (new) |
| conf | 1 ต่ำ / 2 กลาง / 3 สูง | 1 low / 2 medium / 3 high |

## Coordinate systems

Rasters and `.npz` arrays: WGS 84 / UTM zone 47N (EPSG:32647), 30.66 m cells, 261 × 130. Unit polygons and features: EPSG:32647 unless stated; application assets: EPSG:4326.

## Licences and sources

Data: CC BY 4.0. Code: MIT. Third-party inputs remain under their own terms: Copernicus GLO-30 DEM (© DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018, provided under COPERNICUS by the European Union and ESA); ESA WorldCover 2021 v200 (CC BY 4.0); modified Copernicus Sentinel data 2019–2025; GPM IMERG (NASA). National C values: Land Development Department (2020). Soil data: Horpeancharoen et al. (2022).
