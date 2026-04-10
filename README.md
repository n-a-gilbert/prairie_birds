# Limited capacity of landscape features to buffer grassland bird declines

### Data/code DOI:
__________________________________________________________________________________________________________________________________________

## Abstract
Grassland birds are among the most rapidly declining groups of species, largely due to historical habitat losses that have led to extinction debts. Beyond direct habitat manipulations (e.g., prescribed fire), conservation strategies for grassland birds often adopt a landscape-scale perspective, particularly when guiding decisions on reserve acquisition, prioritizing management across existing reserves, and restoring habitat. Additionally, climate is increasingly recognized as a critical consideration for grassland bird conservation, since temperature influences survival and nest success, while precipitation has an outsized effect on primary productivity and habitat structure in grassland ecosystems. Here, we used a 16-year breeding season avian point count dataset from some of the highest-quality patches of remaining grassland habitat in Minnesota, USA, to 1) estimate occurrence trends for 33 species that use grasslands and 2) assess the influence of two traditional habitat variables (amount of grassland habitat, grassland edge density) and two climate factors (temperature precipitation anomaly) on baseline occurrence and trends. Our model reveals declines across bird species that use grasslands; 4 out of 8 grassland-obligate species declined (with ≥95% confidence; 6/8 declined with ≥68% confidence) whereas 6 out of 25 facultative grassland species declined (with ≥95% confidence). The traditional habitat variables were associated with baseline occurrence of relatively few species (e.g., proportion of grassland habitat showed a clear relationship with baseline occurrence of 4 species, including 2 obligates) and buffered declines for only a handful of species, none of them obligates (e.g., amount of grassland habitat showed a buffering effect for Red-winged Blackbird). In contrast, breeding-season climate anomalies showed strong relationships with occurrence of more species (temperature: 7 species; precipitation: 11 species). However, climate effects on trends were less common (e.g., wetter breeding seasons buffered declines of Savannah Sparrows). Considered holistically, our results suggest that grassland bird declines are apparent even in the highest-quality remaining habitats of the upper midwestern US. Habitat factors may offer limited buffering against future grassland bird declines, especially as temperature and precipitation anomalies exert increasingly strong impacts on population trajectories and climate conditions continue to change. 
 $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ <img src="https://github.com/n-a-gilbert/prairie_birds/blob/main/figures/figure_04.png" width="600" />
 
## Repository Directory

### code
 * [01_calculate_grassland_vars.R](./code/01_calculate_grassland_vars.R) Calculate proportion grassland and edge density from contiguous grassland shapefile
 * [01_calculate_weather_anomaly.R](./code/01_calculate_weather_anomaly.R) Calculate weather anomalies from long-term conditions
 * [02_format_data.R](./code/02_format_data.R) Format data for analysis
 * [03_check_linear_quadratic.R](./code/03_check_linear_quadratic.R) Quick check of linear vs. quadratic effects of climate variables
 * [03_fit_models_hpc.R](./code/03_fit_models_hpc.R) Fit models 
 * [04_figure_01.R](./code/04_figure_01.R) Create Figure 1 (study area map)
 * [04_figure_02.R](./code/04_figure_02.R) Create Figure 2 (trend estimates)
 * [04_figure_03.R](./code/04_figure_03.R) Create Figure 3 (effects of ecological covariates on occurrence and trends)
 * [04_figure_04.R](./code/04_figure_04.R) Create Figure 4 (species heatmap showing covariate effects)
 * [04_figure_s01.R](./code/04_figure_s01.R) Create Figure S1
 * [04_figure_s02.R](./code/04_figure_s02.R) Create Figure S2
 * [04_figure_s03.R](./code/04_figure_s03.R) Create Figure S3
 * [04_figure_s04.R](./code/04_figure_s04.R) Create Figure S4
 * [04_table_s01.R](./code/04_table_s01.R) Create Table S01
### data
 * [spatial](./data/spatial) Folder containing spatial data files
   * [SPICE_contig_grass_NoHoles_20260213.shp](./data/spatial/SPICE_contig_grass_NoHoles_20260213.shp) Shapefile for contiguous grassland habitat (.cpg, .dbf, .prj, .sbn, .sbx, and .shx components also present in this folder)
   * [allsites_SPICE_20201228_multipart.shp](./data/spatial/allsites_SPICE_20201228_multipart.shp) Shapefile with reserve boundaries (.cpg, .dbf, .prj, .sbn, .sbx, and .shx components also present in this folder)
   * [ecs_provinces_of_mn_v99a.shp](./data/spatial/ecs_provinces_of_mn_v99a.shp) Shapefile with ecoregions of Minnesota (.cpg, .dbf, .prj, .sbn, .sbx, and .shx components also present in this folder)
 * [MN_prairie_bird_data.xlsx](./data/MN_prairie_bird_data.xlsx) Grassland bird survey data (format received from collaborators); included here since it is used in scripts to calculate weather anomalies and habitat variables.
 * [habitat_vars.csv](./data/habitat_vars.csv) Point-level habitat variables

  | Column name | Type | Description |
  |-------------|------|-------------|
  | point | character | name of survey point nested within reserve |
  | pgrass_250 | double | proportion grassland within 250 m of point |
  | ed_250 | double | edge density (grassland/non-grassland edge) within 250 m of point |
  | pgrass_3000 | double | proportion grassland within 3000 m of point |
  | ed_3000 | double | edge density (grassland/non-grassland edge) within 3000 m of point |

 * [list_of_data_lists3.RData](./data/list_of_data_lists3.RData) List of data lists for modeling. 8 top-level elements (1 for each model). Each element is a list of six, containing model constants (loop control, etc), data, model code, parameters monitored, initial value function, and model label.
 * [mn_prairie_bird_data_clean.csv](./data/mn_prairie_bird_data_clean.csv) Grassland bird survey data. Data dictionary provided below.

  | Column name | Type | Description |
  |-------------|------|-------------|
  | site | character | name of reserve |
  | point | character | name of survey point nested within reserve |
  | easting | double | UTM easting (EPSG code: 26915) |
  | northing | double | UTM northing (EPSG code: 26915) |
  | year | double | year the survey was conducted |
  | visit | double | which of the replicate visits (1, 2, or 3) |
  | date | double | ordinal date (day-of-the-year) of the survey |
  | obs | character | initials of the observer conducting the survey | 
  | start | double | start time (hours after midnight) of survey |
  | sp | character | 4-letter species code |
  | n | double | number of individuals counted during survey |
  | area | double | area (sq. km) of the reserve |
  | ratio | double | perimeter-to-area ratio of the reserve |
  | open | double | proportion open habitat within 250 m of survey point. This was calculated from the 2016 NLCD using Google Earth Engine; open habitats were considered to include the "barren land", "shrub/scrub", "grassland/herbaceous", "sedge/herbaceous", "pasture/hay", and "emergent herbaceous wetlands" classes |
  | anom | double | breeding season (May and June) precipitation anomaly from 1985–2005 average, calculated from Daymet | 

 * [sp_key.csv](./data/sp_key.csv) Key with species identifiers

   | Column name | Type | Meaning |
   |-------------|------|---------|
   | sp | double | species index used in model |
   | code | character | 4-letter species code |
   | common | character | species common name |

 * [weather_anomaly.csv](./data/weather_anomaly.csv) Weather anomalies

   | Column name | Type | Meaning |
   |-------------|------|---------|
   | point | character | name of survey point nested within reserve |
   | year | double | year |
   | panom_br | double | breeding season (May and June) precipitation anomaly |
   | tanom_br | double | breeding season (May and June) temperature anomaly |
   | panom_nb | double | nonbreeding season (Jan-Apr) precipitation anomaly |
   | tanom_nb | double | nonbreeding season (Jan-Apr) temperature anomaly |
   
### figures
 * [figure_01.png](./figures/figure_01.png) Figure 1
 * [figure_02.png](./figures/figure_02.png) Figure 2
 * [figure_03.png](./figures/figure_03.png) Figure 3
 * [figure_04.png](./figures/figure_04.png) Figure 4
 * [figure_s01.png](./figures/figure_s01.png) Figure S1
 * [figure_s02.png](./figures/figure_s02.png) Figure S2
 * [figure_s03.png](./figures/figure_s03.png) Figure S3
 * [figure_s04.png](./figures/figure_s04.png) Figure S4
 
### results
 * [01_pgrass_250_tanom_br_panom_br.RData](./results/01_pgrass_250_tanom_br_panom_br.RData) Model information and model output for model with pgrass (250 m) and breeding-season climate variables. RData with the following objects:
   * constants: constants (loop control, indices) for model
   * data: data used to fit model
   * out: model output with posterior distributions of model parameters
   * settings: MCMC settings used to fit model
   * model_code: Nimble model code 
 * [02_pgrass_250_tanom_nb_panom_nb.RData](./results/02_pgrass_250_tanom_nb_panom_nb.RData) Model information and model output for model with pgrass (250 m) and nonbreeding-season climate variables
 * [03_pgrass_3000_tanom_br_panom_br.RData](./results/03_pgrass_3000_tanom_br_panom_br.RData) Model information and model output for model with pgrass (3000 m) and breeding-season climate variables
 * [04_pgrass_3000_tanom_nb_panom_nb.RData](./results/04_pgrass_3000_tanom_nb_panom_nb.RData) Model information and model output for model with pgrass (3000 m) and nonbreeding-season climate variables
 * [05_ed_250_tanom_br_panom_br.RData](./results/05_ed_250_tanom_br_panom_br.RData) Model information and model output for model with edge density (250 m) and breeding-season climate variables
 * [06_ed_250_tanom_nb_panom_nb.RData](./results/06_ed_250_tanom_nb_panom_nb.RData) Model information and model output for model with edge density (250 m) and nonbreeding-season climate variables
 * [07_ed_3000_tanom_br_panom_br.RData](./results/07_ed_3000_tanom_br_panom_br.RData) Model information and model output for model with edge density (3000 m) and breeding-season climate variables
 * [08_ed_3000_tanom_nb_panom_nb.RData](./results/08_ed_3000_tanom_nb_panom_nb.RData) Model information and model output for model with edge density (3000 m) and nonbreeding-season climate variables
* [prairie_bird_trends_global_tanom_panom_br22026-02-01.RData](./results/prairie_bird_trends_global_tanom_panom_br22026-02-01.RData) Model information and model output for model with original covariates but quadratic effects of breeding-season climate variables to assess linear vs. quadratic effects
