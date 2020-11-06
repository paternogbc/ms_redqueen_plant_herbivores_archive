# folder data/processed

## Files information

This folder conatins 7 files

| N | file     |      description      |
|:----|----------|:-------------:|
| 1 | full_data.csv | The complete processed comparatived dataset used in the main analysis |
| 2 | data_exclude_repro.csv  | The comparative dataset generated after excluding records on reproductive organs |
| 3 | data_full_records | The comparative dataset generated including all records |
| 4 | phylogenetyic_tree_s1.tre | The study phylogeny generated with scenario S1 |
| 5 | phylogenetyic_tree_s2.tre | The study phylogeny generated with scenario S3 |
| 6 | phylogenetyic_trees_300_s2.tre | A multiphylo file with 300 trees generated with scenario S2 |
| 7 | study_species_list_taxonomy | The list of of species  and thier taxonomic classification |


## 1. full_data.csv

Variables description 

| Column    |      description      | unit |
|-----------|:-------------:|-------:|
| order     | The taxonomic order of the species | |
| family    | The taxonomic family of the species | |
| tip_name  | Binomial species name (genus_epithet) matching the phylogenetic tree | |
| maleness  | Flower maleness. Calculated by dividing the dry biomass of androecium by the dry biomass of androeciem + gynoecium ||
| nspe      | The number of insect species associated to each plant species ||
| nfam      | The number of insect families associated to each plant species ||
| ngui      | The number of insect feeding guilds associated to each plant species ||
| shan      | Shannon diversity index of feeding guilds associated to each plant species ||
| height    | The maximum height of each plant species | (m) |
| sla_imp   | The specific leaf area of each plant species | (mm . mg ^-1) |
| nectar    | The categorical level of nectar offer for each plant species |  |
| color     | The category of flower color for each plant species |  |
| shape     | The category of flower shape for each plant species |  |
| polli     | The category of pollinator group for each plant species |  |
| range     | The area of occupancy for each plant species. Calculated as the number of TK 25 grids out of 3000 |  |
