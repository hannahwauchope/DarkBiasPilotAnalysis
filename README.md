# Dark Bias Pilot Analysis
A pilot analysis of LPI and BioTIME data to assess Weighting and Dark Bias according to levels of habitat loss. This analysis was conducted as part of preliminary testing in preparing a grant: "Just how fast are we losing wildlife".

## Rationale and Question
The Living Planet Database and BioTIME are both global scale databases that give annual abundance estimates for wildlife at discrete locations. These abundance estimates are used to estimate change through time for populations, and to create summarised statistics of overall change in global biodiversity.

To accurately represent global biodiversity changes, these datasets would need to be representative of all populations, and especially of the range of changes occurring in populations (increases or declines). We reason that a major driver of whether a population is increasing or declining is the level of habitat loss at its location, and so asked whether the LPD and BioTIME sample representatively from levels of habitat loss across the world. 

Habitat loss, and levels of sampling, vary by continent and biome, and so we further specify our core question to the following:

<p style="text-align:center;">
**How representative are LPD and BioTIME samples of levels of habitat loss experienced in each biome, in each continent?**
</p>

Habitat loss also varies annually, and so we estimate loss across all years for which data is available.

## Map Data
Biome data was obtained from WWF "Terrestrial Ecoregions of the World" (WWF link to download currently down, also available [here](https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/). We simplified these to a coarser scale, for instance, categorising 'Tropical & Subtropical Moist Broadleaf Forests' to simply 'Tropical Forests' etc - see Code Section 1A for all categorisations.

Continent map data was obtained from [TM World Borders](https://koordinates.com/layer/7354-tm-world-borders-03/). (See Code Section 1B)

## Population Data 
Living Planet Database data was accessed [here](https://www.livingplanetindex.org/data_portal), and BioTIME data accessed [here](https://biotime.st-andrews.ac.uk/download.php). Data was filtered to all records with coordinate locations, and after 1950. See Code Section 2

## Habitat Loss Data
Habitat loss data was obtained from the Land Use Harmonisation project, see [Hurtt et al., 2020](https://gmd.copernicus.org/articles/13/5425/2020/). 

We took data from the [LUH2 v2h release - Historic Data - Transitions](https://luh.umd.edu/data.shtml). This gives gridded annual estimates of all transitions between habitat types between 850AD and 2015 (for example 'Primary Forested Land to Secondary Non Forested Land', 'Rangeland to Urban', 'Pasture to Biofuel Cropland'). There are two types of Primary Habitat (Primary Forested Land, Primary Non-forested Land), we considered any transition from these classifications to another as 'Primary Habitat Loss', and summed all transition estimates per grid cell, per year (for all years from 1950-2015). See Code Section 3. 

## Workflow
We extracted Biome classification and continent for each record in LPI and BioTIME (Code Section 4)

Then, for each Biome in each continent we obtained (Code Section 5):
- Estimates of Primary Habitat Loss for every gridcell in every year
- Estimates of Primary Habitat Loss in gridcell-years sampled by LPD and BioTIME

Finally, we created density plots visualising the distribution of these estimates per Biome/Continent grouping (Code Section 6). We could not include all examples in our main figure, see 'Figures' for all Continent x Biome combinations.
