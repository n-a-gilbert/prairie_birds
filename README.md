# Limited capacity of landscape features to buffer grassland bird declines

### Data/code DOI:
__________________________________________________________________________________________________________________________________________

## Abstract
Grassland species are one of the most severely declining groups of birds, largely due to historical habitat loss which has incurred extinction debts. Beyond direct habitat manipulations (e.g., prescribed fire), conservation efforts for grassland birds often take a landscape perspective, particularly when making decisions about acquiring reserves and allocating management efforts among existing reserves. Additionally, climate is increasingly recognized as a critical consideration for grassland bird conservation, since precipitation has a strong effect on habitat structure and grassland habitats generally experience less microclimate buffering compared to forested habitats. Here, we used a 15-year monitoring dataset of breeding bird surveys from the highest-quality patches of remaining grassland habitat in Minnesota, USA, to 1) estimate occurrence trends for 33 species that use grasslands and 2) assess the influence of three traditional landscape factors (reserve size, reserve shape, and amount of open habitat) and one climate factor (precipitation anomaly) on baseline occurrence and trends. Our model revealed declines across the grassland bird community, including for 7 out of 8 grassland obligate species and 3 out of the 25 facultative grassland species. The traditional landscape variables were associated with baseline occurrence of relatively few species (e.g., amount of open habitat showed a clear relationship with baseline occurrence of 3 species, all obligates) and buffered declines for only a handful of species (e.g., amount of open habitat showed a clear buffering effect for only two species, Bobolink and Red-winged Blackbird). In contrast, breeding-season precipitation anomaly showed clear relationships with baseline occurrence for 15 species (with a nearly even mix of positive and negative effects); however, effects on trends were few (wetter conditions buffered declines for Field and Savannah Sparrows but accelerated declines for Dickcissels). Considered holistically, our results indicate that grassland bird declines are apparent even in the highest-quality remaining habitats. Landscape factors may have limited capacity to buffer against grassland bird declines in the future, especially as precipitation anomalies have a larger impact on population trajectories and climate conditions continue to change. 
 $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ <img src="https://github.com/n-a-gilbert/prairie_birds/blob/main/figures/figure_04.png" width="600" />
 
## Repository Directory

### code
 * [figure_01.R](./code/figure_01.R) Create Figure 1 (study area map)
 * [figure_02.R](./code/figure_02.R) Create Figure 2 (trend estimates)
 * [figure_03.R](./code/figure_03.R) Create Figure 3 (effects of ecological covariates on occurrence and trends)
 * [figure_04.R](./code/figure_04.R) Create Figure 4 (species heatmap showing covariate effects)
 * [figure_s01_s03_s05_s08.R](./code/figure_s01_s03_s05_s08.R) Create Figures S01, S03, S05, S08
 * [figure_s02_s04_s06_s09.R](./code/figure_s02_s04_s06_s09.R) Create Figures S02, S04, S06, S09
 * [figure_s07.R](./code/figure_s07.R) Create Figure S07
 * [format_data_fit_model.R](./code/format_data_fit_model.R) Code to format data and fit multispecies occupancy model using Nimble
 * [table_s01.R](./code/table_s01.R) Create Table S01
### data
 * [spatial](./data/spatial) Folder containing spatial data files
   * [allsites_SPICE_20201228_multipart.shp](./data/spatial/allsites_SPICE_20201228_multipart.shp) Shapefile with reserve boundaries (.cpg, .dbf, .prj, .sbn, .sbx, and .shx components also present in this folder)
   * [ecs_provinces_of_mn_v99a.shp](./data/spatial/ecs_provinces_of_mn_v99a.shp) Shapefile with ecoregions of Minnesota (.cpg, .dbf, .prj, .sbn, .sbx, and .shx components also present in this folder)
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
 * 
### figures
### results
