##preparation for jonas paper
source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")
library(reshape2)
# inputDB <- read_rds("N:/COVerAGE-DB/Data/InputDB.rds") %>% 
#   osf_retrieve_file("9dsfk") %>%
#   osf_download(conflicts = "overwrite")
freesz  <- memuse::Sys.meminfo()$freeram@size
n.cores <- 10

setwd('.')
Offsets     <- readRDS("./tmp/harmonized_population.rds") %>%  
  mutate(Age = as.numeric(substr(id, 12,14)),
         Year = substr(id, 8,11),
         Region = substr(id, 1,6),
         Sex = substr(id, 7,7)) %>%   
  group_by(Region, Year, population_source, Age) %>% 
  summarise(population_midyear = sum(population_midyear)) %>% 
  filter(Year == "2020",
         !is.na(population_midyear))
Offsets$Sex <- "b"
Offsets$Region <- gsub("----", "", Offsets$Region)
names(Offsets)[1] <- "Country"
names(Offsets)[5] <- "Population"
Offsets$Region <- "All"

Offsets <- Offsets %>% 
  mutate(Country = case_when(
      Country == "AT" ~ "Austria",
      Country == "BE" ~ "Belgium",
      Country == "BG" ~ "Bulgaria",
      Country == "CL" ~ "Chile",
      Country == "HR" ~ "Croatia",
      Country == "CZ" ~ "Czech Republic",
      Country == "DK" ~ "Denmark",
      Country == "EE" ~ "Estonia",
      Country == "FI" ~ "Finland",
      Country == "FR" ~ "France",
      Country == "DE" ~ "Germany",
      Country == "HU" ~ "Hungary",
      Country == "IS" ~ "Iceland",
      Country == "IT" ~ "Italy",
      Country == "LT" ~ "Lithuania",
      Country == "NL" ~ "Netherlands",
      Country == "NO" ~ "Norway" ,
      Country == "PL" ~ "Poland",
      Country == "PT" ~ "Portugal",
      Country == "SK" ~ "Slovakia",
      Country == "SI" ~ "Slovenia",
      Country == "ES" ~ "Spain",
      Country == "SE" ~ "Sweden",
      Country == "CH" ~ "Switzerland",
      Country == "GB-EAW" ~ "England and Wales",
      Country == "GB-NIR" ~ "Northern Ireland",
      Country == "GB-SCT" ~ "Scotland",
      Country == "US" ~ "USA"
      ))
Offsets <- Offsets[,-c(2,3)]
Offsets$AgeInt <- 1L

#####adding age groups up to 105

Offsets <-
  Offsets %>% 
  mutate(AgeInt = ifelse(AgeInt == 0, 1, AgeInt))


# Split offsets by country/region/sex
oL <-split(Offsets, 
           list(Offsets$Country,Offsets$Region,Offsets$Sex), 
           drop = TRUE)

# Parallelized harmonization
oL1 <- mclapply(
  oL,
  try_step,
  process_function = harmonize_offset_age_p,
  byvars = c("Country","Region","Sex"),
  mc.cores=n.cores )

# Combine offsets in data frame
Offsets <-
  oL1 %>% 
  rbindlist() %>% 
  as.data.frame()


#############Vaccination data from COVerAGE-DB#######################
osf_retrieve_file("9dsfk") %>%
  osf_download(conflicts = "overwrite")

inputDB <-  read_csv("inputDB.zip",
                     skip = 1,
                     col_types = "cccccciccdc") %>% 
  #inputDB <- read.csv("N:/COVerAGE-DB/Data/InputDB.csv", skip = 1) %>% 
  mutate(Date = dmy(Date)) %>% 
  filter(Country %in% c("Czechia",
                        "Northern Ireland", "Iceland", "Sweden", "Austria", "Slovenia", "Germany", "Finland", "Denmark", "Italy",
                        "Bulgaria", "Hungary", "Poland", "Lithuania", "Estonia", "Spain", "Belgium",
                        "France", "Scotland", "Portugal", "Chile", "USA", "Norway", "Switzerland", "Croatia", "Slovakia")) %>% 
  filter(Measure != "Deaths",
         Measure != "Tests",
         Measure != "Cases",
         Short != "US_All") 

#estonia <- subset(inputDB2, Country == "Estonia")

###sumup to total country
czechia <- inputDB %>% 
  filter(Country == "Czechia") %>% 
  group_by(Date, Country, Age, Measure, AgeInt, Metric, Sex) %>% 
  summarise(Value = sum(Value)) %>% 
  mutate(Region = "All",
         Country = "Czech Republic")

###taking the other rds file for countries that did not make it into coverage
 
  scotland <- read_rds("N:/COVerAGE-DB/Automation/Hydra/Scotland_Vaccine.rds") %>% 
    mutate(Date = dmy(Date),
           Value = as.numeric(Value)) 


 usa <- read_rds("N:/COVerAGE-DB/Automation/Hydra/USA_Vaccine.rds") %>% 
   mutate(Date = dmy(Date),
          Value = as.numeric(Value))

netherlands <- read_rds("N:/COVerAGE-DB/Automation/Hydra/Netherlands_Vaccine.rds") %>% 
  mutate(Date = dmy(Date),
         Value = as.numeric(Value))

inputDB3 <- inputDB %>% 
  filter(Country != "Czechia")
inputDB4 <- bind_rows(inputDB3,czechia, scotland, netherlands, usa) %>% 
  # filter(Age != "UNK",
  #        Age != "TOT")%>% 
  filter(Measure != "Vaccinations") %>% 
  filter(Region == "All")


####harmonize metrics

logfile <- ("./tmp/buildlog2.md")

### Get data ########################################################




# this script transforms the inputDB as required, and produces standardized measures and metrics

icolsIN <- colnames(inputDB4)
icols   <- icolsIN[icolsIN != "AgeInt"]

### Remove unnecessary rows #########################################

Z <-
  inputDB4 %>% 
  # TR: This removes sex-totals that are fractions.
  filter(!(Age == "TOT" & Metric == "Fraction"),
         !(Age == "UNK" & Value == 0),
         !(Sex == "UNK" & Value == 0),
         !is.na(Value)) %>% 
  as.data.table() %>% 
  mutate(AgeInt = as.integer(AgeInt))


### Covert UNK Sex, UNK Age to both-sex TOT entries ################

# Log

AA <- Z[ , try_step(
  process_function = resolve_UNKUNK,
  chunk = .SD,
  byvars = c("Country","Region","Date","Measure"),
  logfile = logfile),
  by = list(Country, Region, Date, Measure),
  .SDcols = icols][,..icols]


### Convert fractions ###############################################

# Log 

# Convert sex-specific fractions to counts
A <- AA[ , try_step(process_function = convert_fractions_sexes,
                    chunk = .SD,
                    byvars = c("Code","Measure"),
                    logfile = logfile),
         by = list(Country, Region, Date, Measure), 
         .SDcols = icols][,..icols]

# Convert fractions within sexes to counts
A <- A[ , try_step(process_function = convert_fractions_within_sex,
                   chunk = .SD,
                   byvars = c("Code","Sex","Measure"),
                   logfile = logfile),
        by=list(Country, Region, Date, Sex, Measure), 
        .SDcols = icols][,..icols]

### Distribute counts with unknown age ##############################

# Log

B <- A[ , try_step(process_function = redistribute_unknown_age,
                   chunk = .SD,
                   byvars = c("Code","Sex","Measure"),
                   logfile = logfile),
        by = list(Country, Region, Date, Sex, Measure),
        .SDcols = icols][,..icols]

# ### Scale to totals (within sex) ####################################
#
# # Log
#
C <- B[ , try_step(process_function = rescale_to_total,
                   chunk = .SD,
                   byvars = c("Code","Sex","Measure"),
                   logfile = logfile),
        by = list(Country, Region, Date, Sex, Measure),
        .SDcols = icols][,..icols]

### Derive counts from deaths and CFRs ##############################

# Log
D <- C[ , try_step(process_function = infer_cases_from_deaths_and_ascfr,
                   chunk = .SD,
                   byvars = c("Code", "Sex"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Sex), 
        .SDcols = icols][,..icols]

# Infer deaths from cases and CFRs ##################################

# Log

E <- D[ , try_step(process_function = infer_deaths_from_cases_and_ascfr,
                   chunk = .SD,
                   byvars = c("Code", "Sex"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Sex), 
        .SDcols = icols][,..icols]

# Drop ratio (just to be sure, above call probably did that)
E <- E[Metric != "Ratio"]

### Distribute cases with unkown sex ################################


G <- E[ , try_step(process_function = redistribute_unknown_sex,
                   chunk = .SD,
                   byvars = c("Code", "Age", "Measure"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Age, Measure), 
        .SDcols = icols][,..icols]

### Scale sex-specific data to match combined sex data ##############

# Log

H <- G[ , try_step(process_function = rescale_sexes,
                   chunk = .SD,
                   byvars = c("Code", "Measure"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Measure), 
        .SDcols = icols][,..icols]

# Remove sex totals
H <- H[Age != "TOT"]

### Both sexes combined calculated from sex-specifc #################

# Log

I <- H[ , try_step(process_function = infer_both_sex,
                   chunk = .SD,
                   byvars = c("Code", "Measure"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Measure), 
        .SDcols = icols][,..icols]


### Adjust closeout age #############################################

# Log

J <- I[ , Age := as.integer(Age), ][, ..icols]

# Adjust
J <- J[ , try_step(process_function = maybe_lower_closeout,
                   chunk = .SD, 
                   byvars = c("Code", "Sex", "Measure"),
                   OAnew_min = 85,
                   Amax = 104,
                   logfile = logfile), 
        by = list(Country, Region, Date, Sex, Measure),
        .SDcols = icols][,..icols]


# Formatting 

inputCounts <- J[ , AgeInt := add_AgeInt(Age, omega = 105),
                  by = list(Country, Region, Date, Sex, Measure)][, ..icolsIN] %>% 
  arrange(Country, Region, Sex, Measure, Age) %>% 
  as.data.frame()
### Saving ##########################################################

######Age-harmonizing################################################




inputCounts <- inputCounts %>% 
  arrange(Country, Region, Date, Measure, Sex, Age) %>% 
  group_by(Code, Sex, Measure, Date) %>% 
  mutate(id = cur_group_id(),
         core_id = sample(1:n.cores,size=1,replace = TRUE)) %>% 
  filter(Sex == "b") %>% 
  ungroup() 

# nr rows per core
# inputCounts$core_id %>% table()

ids_in <- inputCounts$id %>% unique() %>% sort()

# Number of subsets per core
# tapply(inputCounts$id,inputCounts$core_id,function(x){x %>% unique() %>% length()})
# Split counts into big chunks
iL <- split(inputCounts,
            inputCounts$core_id,
            drop = TRUE)

### Age harmonization: 5-year age groups ############################


#print(object.size(iL),units = "Mb")
# 5 feb 2021 800 Mb
# length(iL)
tic()
# Apply PCLM to split into 5-year age groups
vaccination_harmonized <- parallelsugar::mclapply(
  iL, 
  harmonize_age_p_bigchunks,
  Offsets = Offsets, # 2.1 Mb data.frame passed to each process
  N = 5,
  OAnew = 100,
  lambda = 1e5,
  mc.cores = n.cores)
toc()


out <- rbindlist(vaccination_harmonized) 
ids_out  <- out$id %>% unique() %>% sort()
failures <- ids_in[!ids_in %in% ids_out]

HarmonizationFailures <-
  inputCounts %>% 
  filter(id %in% failures)




write_rds(out, "./tmp/vaccination_harmonized2.rds")
write_rds(inputCounts, "./tmp/inputCounts2.rds")
write_rds(HarmonizationFailures, "./tmp/failiors.rds")




harmonized <- out
failed <- HarmonizationFailures
failed <- failed[,-12]



over60 <- harmonized %>% 
  filter(Age >= 60,
         Sex == "b") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value))

over60.2 <- failed %>% 
  filter(Age >= 60,
         Sex == "b",
         Country != "Netherlands") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value))

over60.3 <- rbind(over60, over60.2)

over60.pop <- Offsets%>% 
  filter(Age >= 60,
         Sex == "b",
         Region == "All") %>% 
  group_by(Country) %>% 
  summarise(Population = sum(Population))


over60.rates <- merge(over60.3, over60.pop)
over60.rates$rate <- over60.rates$Value / over60.rates$Population

write_rds(over60.rates, "./tmp/rates_v2.rds")

totals <- harmonized %>% 
  filter(Sex == "b") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value))

totals.2 <- failed %>% 
  filter(Sex == "b") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value)) %>% 
  filter(Country != "Portugal",
         Country != "Norway")

totals.3 <- rbind(totals, totals.2)

totals.pop <- Offsets%>% 
  filter(Sex == "b",
         Region == "All") %>% 
  group_by(Country) %>% 
  summarise(Population = sum(Population))


totals.rates <- merge(totals.3, totals.pop)
totals.rates$rate <- totals.rates$Value / totals.rates$Population

write_rds(totals.rates, "./tmp/totals_rates_v2.rds")


under60 <- harmonized %>% 
  filter(Age < 60,
         Sex == "b") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value))

under60.2 <- failed %>% 
  filter(Age < 60,
         Sex == "b",
         Country != "Netherlands") %>% 
  group_by(Date, Country, Region, Measure) %>% 
  summarise(Value = sum(Value))

under60.3 <- rbind(under60, under60.2)

under60.pop <- Offsets%>% 
  filter(Age < 60,
         Sex == "b",
         Region == "All") %>% 
  group_by(Country) %>% 
  summarise(Population = sum(Population))


under60.rates <- merge(under60.3, under60.pop)
under60.rates$rate <- under60.rates$Value / under60.rates$Population

write_rds(under60.rates, "./tmp/rates_v2_under60.rds")

