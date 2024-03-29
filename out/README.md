Essential output files of the analysis are archived with Zenodo. See https://doi.org/10.5281/zenodo.6241024 for the files

- 40-lifetables.csv
- 30-lt_input.csv

from which all the results in the paper can be derived.

## 40-lifetables.csv

Life table statistics 2015 through 2021 by sex, region and quarter with uncertainty quantiles based on Poisson replication of death counts. Actual life tables and expected life tables (under the assumption of pre-COVID mortality trend continuation) are provided.

## 30-lt_input.csv

Life table input data.

- `id`: unique row identifier
- `region_iso`: iso3166-2 region codes
- `sex`: Male, Female, Total
- `year`: iso year
- `age_start`: start of age group
- `age_width`: width of age group, Inf for age_start 100, otherwise 1
- `nweeks_year`: number of weeks in that year, 52 or 53
- `death_total`: number of deaths by any cause
- `population_py`: person-years of exposure (adjusted for leap-weeks and missing weeks in input data on all cause deaths)
- `death_total_nweeksmiss`: number of weeks in the raw input data with at least one missing death count for this region-sex-year stratum. missings are counted when the week is implicitly missing from the input data or if any NAs are encounted in this week or if age groups are implicitly missing for this week in the input data (e.g. 40-45, 50-55)
- `death_total_minnageraw`: the minimum number of age-groups in the raw input data within this region-sex-year stratum
- `death_total_maxnageraw`: the maximum number of age-groups in the raw input data within this region-sex-year stratum
- `death_total_minopenageraw`: the minimum age at the start of the open age group in the raw input data within this region-sex-year stratum
- `death_total_maxopenageraw`: the maximum age at the start of the open age group in the raw input data within this region-sex-year stratum
- `death_total_source`: source of the all-cause death data
- `death_total_prop_q1`: observed proportion of deaths in first quarter of year
- `death_total_prop_q2`: observed proportion of deaths in second quarter of year
- `death_total_prop_q3`: observed proportion of deaths in third quarter of year
- `death_total_prop_q4`: observed proportion of deaths in fourth quarter of year
- `death_expected_prop_q1`: expected proportion of deaths in first quarter of year
- `death_expected_prop_q2`: expected proportion of deaths in second quarter of year
- `death_expected_prop_q3`: expected proportion of deaths in third quarter of year
- `death_expected_prop_q4`: expected proportion of deaths in fourth quarter of year
- `population_midyear`: midyear population (July 1st)
- `population_source`: source of the population count/exposure data
- `death_covid`: number of deaths due to covid
- `death_covid_date`: number of deaths due to covid as of <date>
- `death_covid_nageraw`: the number of age groups in the covid input data
- `ex_wpp_estimate`: life expectancy estimates from the World Population prospects for a five year period, merged at the midpoint year
- `ex_hmd_estimate`: life expectancy estimates from the Human Mortality Database
- `nmx_hmd_estimate`: death rate estimates from the Human Mortality Database
- `nmx_cntfc`: Lee-Carter death rate projections based on trend in the years 2015 through 2019

### Deaths

- source:
  - STMF input data series (https://www.mortality.org/Public/STMF/Outputs/stmf.csv)
  - ONS for GB-EAW pre 2020
  - CDC for US pre 2020
- STMF:
  - harmonized to single ages via pclm
    - pclm iterates over country, sex, year, and within-year age grouping pattern and converts irregular age groupings, which may vary by country, year and week into a regular age grouping of 0:110
    - smoothing parameters estimated via BIC grid search seperately for every pclm iteration
    - last age group set to [110,111)
    - ages 100:110+ are then summed into 100+ to be consistent with mid-year population information
  - deaths in unknown weeks are considered; deaths in unknown ages are not considered
- ONS:
  - data already in single ages
  - ages 100:105+ are summed into 100+ to be consistent with mid-year population information
  - PCLM smoothing applied to for consistency reasons
- CDC:
  - The CDC data comes in single ages 0:100 for the US. For 2020 we only have the STMF data in a much coarser age grouping, i.e. (0, 1, 5, 15, 25, 35, 45, 55, 65, 75, 85+). In order to calculate life-tables in a manner consistent with 2020, we summarise the pre 2020 US death counts into the 2020 age grouping and then apply the pclm ungrouping into single year ages, mirroring the approach to the 2020 data

### Population

- source:
  - for years 2000 to 2019: World Population Prospects 2019 single year-age population estimates 1950-2019
  - for year 2020: World Population Prospects 2019 single year-age population projections 2020-2100
- mid-year population
- mid-year population translated into exposures:
  - if a region reports annual deaths using the Gregorian calendar
    definition of a year (365 or 366 days long) set exposures equal
    to mid year population estimates
  - if a region reports annual deaths using the iso-week-year
    definition of a year (364 or 371 days long), and if there is a
    leap-week in that year, set exposures equal to
    371/364\*mid_year_population to account for the longer reporting
    period. in years without leap-weeks set exposures equal
    to mid year population estimates. further multiply by fraction of
    observed weeks on all weeks in a year.

### COVID deaths

- source: COVerAGE-DB (https://osf.io/mpwjq/)
- the data base reports cumulative numbers of COVID deaths over days of a year, we extract the most up to date yearly total

### External life expectancy estimates

- source:
  - World Population Prospects (https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_Life_Table_Medium.csv), estimates for the five year period 2015-2019
  - Human Mortality Database (https://mortality.org/), single year and age tables
