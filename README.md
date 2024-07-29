# Differential responses to weather and habitat conditions explain spatial variation in winter abundance trends in a migratory bird of conservation concern. *Journal of Applied Ecology*.
## Schindler, A. R., A. D. Fox, C. K. Wikle, B. M. Ballard, A. J. Walsh, S. B. A. Kelly, M. D. Weegman. 

#### Model code and associated data files to quantify geographic patterns of abundance trends among Greenland white-fronted geese (*Anser albifrons flavirostris*) wintering sites and examine how weather and land-cover changes explain variation in wintering site abundance. 
___
### Authors
Alexander R. Schindler  
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada

Anthony D. Fox  
Department of Ecoscience, Aarhus University, Aarhus, Denmark

Christopher K. Wikle  
Department of Statistics, University of Missouri, Columbia, MO, USA

Bart M. Ballard  
Caesar Kleberg Wildlife Research Institute, Texas A&M University-Kingsville, Kingsville, TX, USA

Alyn J. Walsh  
National Parks and Wildlife Service, Dublin, Ireland

Seán B. A. Kelly  
National Parks and Wildlife Service, Dublin, Ireland

Mitch D. Weegman  
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada
___
### Code files
- `01_estimating_site_trends.R`: code for estimating abundance trends of individual Greenland white-fronted goose wintering sites
- `02_determining_number_of_groups.R` : code for exploratory analysis to determine number of groups to use in the latent class analysis
- `03_LCA.R`: code for the latent class analysis to group wintering sites into groups based on abundance trends and location (i.e., latitude and longitude) 
- `04_fix_label_switching` : code to fix label-switching in LCA model results
- `05_ss_env.R`: code for estimating the effects of environmental variables on annual trends in abundance

### Data files
- `count_data.csv`: Greenland white-fronted goose wintering site count data
- `roost_site_coordinates.csv`: (scaled) coordinates of roosting locations for all Greenland white-fronted goose wintering sites
- `pland_data.csv`: percent land cover data for all Greenland white-fronted goose wintering sites
- `weather_data.csv`: weather data 
- `winter_GDD_data.csv` : cumulative growing degree days (GDD) for all Greenland white-fronted goose wintering sites

### Data file column names (see manuscript for further details)
- `site_name`: name of each Greenland white-fronted goose wintering site
- `site_id`: numerical identifier for each Greenland white-fronted goose wintering site
- `year`: year in the study (1984-2018)
- `year_id`: numerical identifier for each year in the study (i.e., 1-35 corresponds to 1984-2018)
- `count`: number of Greenland white-fronted geese counted
- `lat_scale`: (scaled) latitude of the corresponding wintering flock
- `lon_scale`: (scaled) longitude of the corresponding wintering flock
- `composite_year` : 5-year composite range of the corresponding land cover data
- `grass`: % grassland cover within a 15km radius of the corresponding wintering flock
- `cereal`: % cereal crop cover within a 15km radius of the corresponding wintering flock
- `peat_bog`: % peat bog cover within a 15km radius of the corresponding wintering flock
- `scaled_grass`: (scaled) % grassland cover within a 15km radius of the corresponding wintering flock
- `scaled_cereal`: (scaled) % cereal crop cover within a 15km radius of the corresponding wintering flock
- `scaled_peat_bog`: (scaled) % peat bog cover within a 15km radius of the corresponding wintering flock
- `spring_storm_days`: number of days with a severe storm during spring migration
- `spring_precip`: cumulative precipitation during spring staging
- `pre_breed_freeze`: number of days below freezing during the pre-breeding period
- `post_breed_precip`: cumulative precipitation during the post breeding period
- `autumn_storm_days`: number of days with a severe storm during autumn migration
- `autumn_freeze`: number of days below freezing during autumn staging
- `scaled_spring_storm_days`: (scaled) number of days with a severe storm during spring migration
- `scaled_spring_precip`: (scaled) cumulative precipitation during spring staging
- `scaled_pre_breed_freeze`: (scaled) number of days below freezing during the pre-breeding period
- `scaled_post_breed_precip`: (scaled) cumulative precipitation during the post breeding period
- `scaled_autumn_storm_days`: (scaled) number of days with a severe storm during autumn migration
- `scaled_autumn_freeze`: (scaled) number of days below freezing during autumn staging
- `winter_GDD` : cumulative growing degree days (GDD) on wintering areas from 01 January to day of goose departure from wintering areas
- `scaled_winter_GDD` : (scaled) cumulative growing degree days (GDD) on wintering areas from 01 January to day of goose departure from wintering areas
