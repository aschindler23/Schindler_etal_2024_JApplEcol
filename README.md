# Habitat conditions during winter explain metapopulation dynamics and a range shift in a long-distance migratory bird. *Journal of Applied Ecology*.
## Schindler, A. R., A. D. Fox, C. K. Wikle, B. M. Ballard, A. J. Walsh, S. B. A. Kelly, M. D. Weegman. 

#### Model code and associated data files to identify common subpopulation abundance trends and associated environmental drivers of these trends in Greenland white-fronted geese (*Anser albifrons flavirostris*). 
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
### Files
- `estimate_flock_trends.R`: code for estimating mean trends in flock abundance  
- `latent_class_analysis.R`: code for the latent class analysis to determine common trends among flocks  
- `estimate_env_effects_on_lambda.R`: code for estimating the effects of environmental variables on annual trends in abundance
- `flock_counts.csv`: flock-level Greenland white-fronted goose count data
- `flock_coordinates.csv`: coordinates of Greenland white-fronted goose wintering flocks
- `land_cover_data.csv`: wintering land cover data
- `weather_data.csv`: weather data 

### Data file column names (see manuscript for further details)
- `flock`: numerical identifier for each Greenland white-fronted goose wintering flock
- `year`: numerical identifier for each year in the study (i.e., 1-35 corresponds to 1984-2018)
- `count`: number of Greenland white-fronted geese counted
- `lat_scale`: (scaled) latitude of the corresponding wintering flock
- `lon_scale`: (scaled) longitude of the corresponding wintering flock
- `grass`: % grassland cover within a 15km radius of the corresponding wintering flock
- `cereal`: % cereal crop cover within a 15km radius of the corresponding wintering flock
- `peat_bog`: % peat bog cover within a 15km radius of the corresponding wintering flock
- `ndvi`: % average NDVI within a 15km radius of the corresponding wintering flock
- `spring_storm_days`: number of days with a severe storm during spring migration
- `spring_precip`: cumulative precipitation during spring staging
- `pre_breed_freeze`: number of days below freezing during the pre-breeding period
- `post_breed_precip`: cumulative precipitation during the post breeding period
- `fall_storm_days`: number of days with a severe storm during fall migration
- `fall_freeze`: number of days below freezing during fall staging

