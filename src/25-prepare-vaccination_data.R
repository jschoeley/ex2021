##preparation of vaccination data
source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")
library(reshape2)
library(tidyverse)
library(readr)
library(lubridate)
library(ggplot2)
library(osfr)
library(covidAgeData)


freesz  <- memuse::Sys.meminfo()$freeram@size
n.cores <- 10

setwd('.')

##getting offsets used by COVerAGE DB#################
Offsets <- read.csv("https://files.de-1.osf.io/v1/resources/mpwjq/providers/osfstorage/5ef371ed65982802b2cf2220?action=download&direct&version=5", skip=1, header=TRUE)

#############Vaccination data from COVerAGE-DB#######################
#osf_retrieve_file("9dsfk")
# osf_download(conflicts = "overwrite")

#https://files.de-1.osf.io/v1/resources/mpwjq/providers/osfstorage/5f3ed659746a8100ad1a2420?action=download&direct&version=387
cdb_repo <- osf_retrieve_node("mpwjq")
osf_ls_files(cdb_repo, path = "Data") %>%
  dplyr::filter(name == "inputDB.zip") %>%
  osf_download(path= "./dat/coverage/",conflicts = "overwrite")


inputDB <-  read_csv("./dat/coverage/inputDB.zip",
                     skip = 1,
                     col_types = "cccccciccdc") %>% 
  #inputDB <- read.csv("N:/COVerAGE-DB/Data/InputDB.csv", skip = 1) %>% 
  mutate(Date = dmy(Date)) %>% 
  filter(Country %in% c("Czechia", "Northern Ireland", "Iceland", "Sweden", "Austria", "Slovenia", "Germany", "Finland", "Denmark", "Italy",
                        "Bulgaria", "Hungary", "Poland", "Lithuania", "Estonia", "Spain", "Belgium","France", "Scotland", "Portugal", "Chile", 
                        "USA", "Norway", "Switzerland", "Croatia", "Slovakia", "Netherlands")) %>% 
  filter(Measure != "Deaths",
         Measure != "Tests",
         Measure != "Cases")%>% 
  mutate(Country = case_when(
    Country == "Czechia" ~ "Czech Republic",
    TRUE ~ Country
  )) %>% 
  filter(Measure != "Vaccinations") %>% 
  filter(Region == "All")

##Lithuania needs special handling, so the harmonization can work
#Lit <- inputDB %>% 
#  filter(Country == "Lithuania")

####harmonize metrics

logfile <- ("./tmp/buildlog2.md")


# this script transforms the data as required, and produces standardized measures and metrics

icolsIN <- colnames(inputDB)
icols   <- icolsIN[icolsIN != "AgeInt"]

### Remove unnecessary rows #########################################

Z <-
  inputDB %>% 
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



B <- A[ , try_step(process_function = redistribute_unknown_age,
                   chunk = .SD,
                   byvars = c("Code","Sex","Measure"),
                   logfile = logfile),
        by = list(Country, Region, Date, Sex, Measure),
        .SDcols = icols][,..icols]

# ### Scale to totals (within sex) ####################################

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


I <- H[ , try_step(process_function = infer_both_sex,
                   chunk = .SD,
                   byvars = c("Code", "Measure"),
                   logfile = logfile), 
        by = list(Country, Region, Date, Measure), 
        .SDcols = icols][,..icols]


### Adjust closeout age #############################################



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

######Age-harmonizing################################################

inputCounts <- inputCounts %>% 
  arrange(Country, Region, Date, Measure, Sex, Age) %>% 
  group_by(Code, Sex, Measure, Date) %>% 
  mutate(id = cur_group_id(),
         core_id = sample(1:n.cores,size=1,replace = TRUE)) %>% 
  filter(Sex == "b") %>% 
  ungroup() 


ids_in <- inputCounts$id %>% unique() %>% sort()

# Split counts into big chunks
iL <- split(inputCounts,
            inputCounts$core_id,
            drop = TRUE)

### Age harmonization: 5-year age groups ############################

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

###################Calculation rates for total population, 60+ and under 60######

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

