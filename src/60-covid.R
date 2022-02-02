source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")
library(dplyr)
library(ggrepel)

##reading in vaccination data
setwd('.')

rates <- read_rds("./dat/coverage/rates_v2.rds") %>%  
  filter(Country != "Northern Ireland")
###fill in missing dates and select dates and vaccination status

vacc2 <- rates %>% 
  filter(Measure == "Vaccination2") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`) %>% 
  ungroup() %>% 
  tidyr::complete(Date, nesting(Country), fill=list(rate=0))


vacc3 <- rates %>% 
  filter(Measure == "Vaccination3") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`)

vacc3_nov <- vacc3 %>% 
  filter(Date == "2021-11-30") %>% 
  select(Country, rate)
names(vacc3_nov)[2] <- "rate_v3_nov"



vacc2_august <- vacc2 %>% 
  filter(Date == "2021-08-31") %>% 
  select(Country, rate)
names(vacc2_august)[2] <- "rate_v2_aug"



vacc2_integral <- vacc2 %>%
  filter(Date >= "2021-01-01" , Date <= "2021-12-31" ) %>% 
  group_by(Country) %>% 
  summarise(
    vac_int = rate %>% mean(na.rm = T)) %>% 
  filter(Country != "Northern Ireland",
         Country != "Finland")
names(vacc2_integral)[2] <- "rate_v2_int"

##all vaccine measures
vacc <- merge(vacc2_august, vacc2_integral, all=TRUE)
vacc <- merge(vacc, vacc3_nov, all=TRUE)

# vacc2_integral2 <- vacc2 %>%
#   filter(Date >= "2021-01-01" , Date <= "2021-12-31" ) %>% 
#   arrange(Country, Date) %>% 
#   group_by(Country) %>% 
#   mutate(areabelow = lag(rate, default = 0) + (rate - lag(rate, default = 0))/2) 
#   #mutate(areabelow = cumsum(areabelow)) 
#   summarise(areabelow = sum(areabelow)) %>% 
#   mutate(areabelow = areabelow / 356)



#plotting
data11 <- merge(data, vacc2_integral_65)
data11$ex_diff_q0.5_2021_month <- data11$ex_diff_q0.5_2021 * 12
data11$diff_all <- (data11$ex_diff_q0.5_2020 + data11$ex_diff_q0.5_2021)*12





###prepare the data for ex
data <- as_tibble(fig$e0diff$data)
names(data)[25] <- "Country"

data_new <- read_rds("U:/gits/ex2021/out/e0avgdiff.rds") %>% 
  filter(age == "0",
         sex == "T")


e60 <- as_tibble(dat$lifetables) %>% 
  filter(age == "60", 
         year == "2021",
         sex == "T") %>% 
  select(region_iso, ex_diff_q0.5)
names(e60)[2] <- "e60_diff_q0.5_2021"

e02019 <- as_tibble(dat$lifetables) %>% 
  filter(age == "0", 
         year == "2019",
         sex == "T") %>% 
   select(region_iso, ex_q0.5)


e02021 <- as_tibble(dat$lifetables) %>% 
  filter(age == "0", 
         year == "2021",
         sex == "T") %>% 
  select(region_iso, ex_q0.5)
names(e02021)[2] <- "e0_2021"

# names(e60)[2] <- "e60_diff_q0.5_2021"


# dekompo20_21 <- as_tibble(fig$arriaga$data) %>% 
#   filter(year== "2021",
#          sex == "T",
#          age_group == "[65,75)"|
#            age_group == "[75,85)"|
#            age_group == "[85,Inf)") %>% 
#   group_by(region_iso) %>% 
#   summarise(dekompo_65_20_21 = sum(e0_arriaga_total_q0.5 * 12))
# 
# dekompo19_21 <- as_tibble(fig$arriaga$data) %>% 
#   filter(year== "2021"|
#            year == "2020",
#          sex == "T",
#          age_group == "[65,75)"|
#            age_group == "[75,85)"|
#            age_group == "[85,Inf)") %>% 
#   group_by(region_iso) %>% 
#   summarise(dekompo_65_10_21 = sum(e0_arriaga_total_q0.5) * 12)

data <- merge(data, e60)
#data <- merge(data, dekompo20_21)
#data <- merge(data, dekompo19_21)
data <- merge(data, data_new)
data <- merge(data, e02019)
data <- merge(data, e02021)

data10 <- merge(data, vacc2_integral)




data10 %>% 
  ggplot() +
  geom_point(mapping = aes(x = e0201 , y = e0_2021, size = vac_int))+
  geom_text_repel(mapping = aes(x = e0201 , y = e0_2021, label = region_iso)) +
 # scale_shape_manual(values=1:36) +
  scale_size_area()+
  coord_fixed()+
  geom_abline()+
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake 60+")


data10 %>% 
  ggplot() +
  geom_point(mapping = aes(x = vac_int , y = e0_2021 - e0201))+
  geom_text_repel(mapping = aes(x = vac_int , y =  e0_2021 - e0201, label = region_iso)) +
  # scale_shape_manual(values=1:36) +
  scale_size_area()+
  geom_abline()+
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake 60+")


data$e0201 <- data$ex_q0.5 + (2* data$mean)
###### full vaccination 60+ at the end of august 2021
data2 <- merge(data, vacc2_august)
data2$diff_all <- (data2$ex_diff_q0.5_2020 + data2$ex_diff_q0.5_2021)*12
data2$ex_diff_q0.5_2021_month <- data2$ex_diff_q0.5_2021 * 12
data2$e60_diff_q0.5_2021_month <- data2$e60_diff_q0.5_2021 * 12


##plotting with linear fit line
# ggplot(data2, aes(rate, diff_all, color = Country)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")

data2 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60.png", plot = last_plot(), dpi = 100, width = 7, height = 6)

data2 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = diff_all, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Change in e0 from 2019 to 2021 in months", 
       title = "Change in e0 from 2019 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60_2019.png", plot = last_plot(), dpi = 100, width = 7, height = 6)


##regression modell


model_60_aug_20_21 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data2)
summary(model_60_aug_20_21)

model_60_aug_19_21 <- lm(diff_all ~ rate, data=data2)
summary(model_60_aug_19_21)

###correlations
cor.test(data2$ex_diff_q0.5_2021_month, data2$rate, method = "pearson")
cor.test(data2$ex_diff_q0.5_2021_month, data2$rate, method = "spearman")


cor.test(data2$diff_all, data2$rate, method = "pearson")
cor.test(data2$diff_all, data2$rate, method = "spearman")

model_60_aug_20_21_e60 <- lm(e60_diff_q0.5_2021_month ~ rate, data=data2)
summary(model_60_aug_20_21_e60)
cor.test(data2$e60_diff_q0.5_2021_month, data2$rate, method = "pearson")
cor.test(data2$e60_diff_q0.5_2021_month, data2$rate, method = "spearman")


model_60_aug_20_21_dekompo <- lm(dekompo_65 ~ rate, data=data2)
summary(model_60_aug_20_21_dekompo)
cor.test(data2$dekompo_65, data2$rate, method = "pearson")
cor.test(data2$dekompo_65, data2$rate, method = "spearman")

##ranked data
data2 <- data2 %>% 
  arrange(ex_diff_q0.5_2021) %>% 
  mutate(rank_e0 = row_number())
data2 <- data2 %>% 
  arrange(rate) %>% 
  mutate(rank_vacc = row_number())

##plot ranked data
data2 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rank_vacc, y = rank_e0, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Rank of Percentage Vaccinated 60+ at the end of August", y="Rank of Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/rank_vacc_august_60_2020.png", plot = last_plot(), dpi = 100, width = 7, height = 6)



###dekompositional effekt and vaccination uptake for 65+ at the end of august

data65 <- merge(data, vacc2_65_august)
data65$diff_all <- (data65$ex_diff_q0.5_2020 + data65$ex_diff_q0.5_2021)*12
data65$ex_diff_q0.5_2021_month <- data65$ex_diff_q0.5_2021 * 12
data65$e60_diff_q0.5_2021_month <- data65$e60_diff_q0.5_2021 * 12


data65 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = dekompo_65_20_21, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 65+ at the end of August", y="Decompositional Age Effekt from 2020 to 2021 65+", 
       title = "Decompositional Age Effekt 2020 to 2021 and Vaccination Uptake 65+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60_decompo.png", plot = last_plot(), dpi = 100, width = 7, height = 6)

data65 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = dekompo_65_10_21, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 65+ at the end of August", y="Decompositional Age Effekt from 2019 to 2021 65+", 
       title = "Decompositional Age Effekt 2019 to 2021 and Vaccination Uptake 65+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60_decompo_2019.png", plot = last_plot(), dpi = 100, width = 7, height = 6)



model_dekompo1 <- lm(dekompo_65_20_21 ~ rate, data=data65)
summary(model_dekompo1)

model_dekompo2 <- lm(dekompo_65_10_21 ~ rate, data=data65)
summary(model_dekompo2)

###correlations
cor.test(data65$dekompo_65_20_21, data65$rate, method = "pearson")
cor.test(data65$dekompo_65_20_21, data65$rate, method = "spearman")


cor.test(data65$dekompo_65_10_21, data65$rate, method = "pearson")
cor.test(data65$dekompo_65_10_21, data65$rate, method = "spearman")




####3rd vaccination at the end of november
data9 <- merge(data, vacc3_nov)
data9$ex_diff_q0.5_2021_month <- data9$ex_diff_q0.5_2021 * 12
data9$diff_all <- (data9$ex_diff_q0.5_2020 + data9$ex_diff_q0.5_2021)*12


data9 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Boosterd 60+ at the end of November", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/vacc3_november.png", plot = last_plot(), dpi = 100, width = 7, height = 6)


data9 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = diff_all, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Boosterd 60+ at the end of November", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2019 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/vacc3_november_2019.png", plot = last_plot(), dpi = 100, width = 7, height = 6)


##regression modell
model_60_nov_20_21_3rd <- lm(ex_diff_q0.5_2021_month ~ rate, data=data9)
summary(model_60_nov_20_21_3rd)


###correlations
cor.test(data9$ex_diff_q0.5_2021_month, data9$rate, method = "pearson")
cor.test(data9$ex_diff_q0.5_2021_month, data9$rate, method = "spearman")







#plotting
data10 <- merge(data, vacc2_integral)
data10$ex_diff_q0.5_2021_month <- data10$ex_diff_q0.5_2021 * 12
data10$diff_all <- (data10$ex_diff_q0.5_2020 + data10$ex_diff_q0.5_2021)*12

data10 %>% 
  ggplot() +
  geom_point(mapping = aes(x = vac_int, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:23) +
  labs(x="Integral Vaccination Measure 60+ [0, 1]", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Integral Vaccination Measure")
ggsave("U:/gits/ex2021/tmp/vacc2_integral.png", plot = last_plot(), dpi = 100, width = 7, height = 6)


data10 %>% 
  ggplot() +
  geom_point(mapping = aes(x = vac_int, y = diff_all, color = Country, shape = Country)) +
  scale_shape_manual(values=1:23) +
  labs(x="Integral Vaccination Measure 60+ [0, 1]", y="Change in e0 from 2019 to 2021 in months", 
       title = "Change in e0 from 2019 to 2021 and Integral Vaccination Measure")
ggsave("U:/gits/ex2021/tmp/vacc2_integral_2019.png", plot = last_plot(), dpi = 100, width = 7, height = 6)
  
  

  

