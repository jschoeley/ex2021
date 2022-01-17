library(dplyr)
rates <- read_rds("U:/gits/ex2021/dat/coverage/rates.rds") %>% 
  filter(Country != "Lithuania")
totals.rates <- read_rds("U:/Jonas/totals_rates.rds")%>% 
  filter(Country != "Lithuania")
###fill in missing dates and select dates and vaccination status

vacc2 <- rates %>% 
  filter(Measure == "Vaccination2") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`)

vacc3 <- rates %>% 
  filter(Measure == "Vaccination3") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`)

vacc3_nov <- vacc3 %>% 
  filter(Date == "2021-11-30") %>% 
  select(Country, rate)


totals_vacc2 <- totals.rates %>% 
  filter(Measure == "Vaccination2") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`)

vacc2_august <- vacc2 %>% 
  filter(Date == "2021-08-31") %>% 
  select(Country, rate)

vacc2_september <- vacc2 %>% 
  filter(Date == "2021-09-30") %>% 
  select(Country, rate)

totals_vacc2_september <- totals_vacc2 %>% 
  filter(Date == "2021-09-30") %>% 
  select(Country, rate)

totals_vacc2_august <- totals_vacc2 %>% 
  filter(Date == "2021-08-30") %>% 
  select(Country, rate)

# min30 <- vacc2 %>% 
#   filter(rate >=0.3) %>% 
#   group_by(Country) %>% 
# mutate(ticker = row_number()) %>% 
#   filter(ticker == "1",
#          Country != "Finland",
#          Country != "Northern Ireland")
# 
# min50 <- vacc2 %>% 
#   filter(rate >=0.5) %>% 
#   group_by(Country) %>% 
#   mutate(ticker = row_number()) %>% 
#   filter(ticker == "1",
#          Country != "Finland",
#          Country != "Northern Ireland")


#####only those countries with lost of e0 from 2020 to 2021
vacc2_august_lost <- vacc2 %>% 
  filter(Date == "2021-08-31") %>% 
  select(Country, rate) %>% 
  filter(Country == "Bulgaria"|
         Country == "Slovakia"|
           Country == "USA"|
           Country == "Poland"|
           Country == "Estonia"|
           Country == "Chile"|
           Country == "Hungary"|
           Country == "Czech Republic"|
           Country == "Croatia"|
           Country == "Scotland")



###prepare the data for ex
data <- as_tibble(fig$e0diff$data)
names(data)[25] <- "Country"

e60 <- as_tibble(dat$lifetables) %>% 
  filter(age == "60", 
         year == "2021",
         sex == "T") %>% 
  select(region_iso, ex_diff_q0.5)
names(e60)[2] <- "e60_diff_q0.5_2021"


dekompo <- as_tibble(fig$arriaga$data) %>% 
  filter(year== "2021",
         sex == "T",
         age_group == "[65,75)"|
           age_group == "[75,85)"|
           age_group == "[85,Inf)") %>% 
  group_by(region_iso) %>% 
  summarise(dekompo_65 = sum(e0_arriaga_total_q0.5))

data <- merge(data, e60)
data <- merge(data, dekompo)


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
ggsave("U:/gits/ex2021/tmp/vacc_august_60.png", plot = last_plot(), dpi = 100)

data2 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = diff_all, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Change in e0 from 2019 to 2021 in months", 
       title = "Change in e0 from 2019 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60_2019.png", plot = last_plot(), dpi = 100)

data2 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = dekompo_65, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 60+ at the end of August", y="Dekompositional Effekt 65+", 
       title = "Dekompositional Effekt 2020 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_august_60_dekompo.png", plot = last_plot(), dpi = 100)

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
ggsave("U:/gits/ex2021/tmp/rank_vacc_august_60_2020.png", plot = last_plot(), dpi = 100)







####fully vaccination 60+ at the end of september
data3 <- merge(data, vacc2_september)
data3$ex_diff_q0.5_2021_month <- data3$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data3, aes(rate, ex_diff_q0.5_2021_month, color = Country)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data3 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated 60+ at the end of September", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_september_60.png", plot = last_plot(), dpi = 100)


##regression modell
model_60_sep_20_21 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data3)
summary(model_60_sep_20_21)


###correlations
cor.test(data3$ex_diff_q0.5_2021_month, data3$rate, method = "pearson")
cor.test(data3$ex_diff_q0.5_2021_month, data3$rate, method = "spearman")



####total vaccination at the end of august
data4 <- merge(data, totals_vacc2_august)
data4$ex_diff_q0.5_2021_month <- data4$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data4, aes(rate, ex_diff_q0.5_2021_month, color = country_name)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data4 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated at the end of August", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/total_vacc_august.png", plot = last_plot(), dpi = 100)


##regression modell
model_tot_aug_20_21 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data4)
summary(model_tot_aug_20_21)


###correlations
cor.test(data4$ex_diff_q0.5_2021_month, data4$rate, method = "pearson")
cor.test(data4$ex_diff_q0.5_2021_month, data4$rate, method = "spearman")


####total vaccination at the end of september
data7 <- merge(data, totals_vacc2_september)
data7$ex_diff_q0.5_2021_month <- data7$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data7, aes(rate, ex_diff_q0.5_2021_month, color = country_name)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data7 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated at the end of September", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/total_vacc_september.png", plot = last_plot(), dpi = 100)


##regression modell
model_tot_sep_20_21 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data7)
summary(model_tot_sep_20_21)


###correlations
cor.test(data7$ex_diff_q0.5_2021_month, data7$rate, method = "pearson")
cor.test(data7$ex_diff_q0.5_2021_month, data7$rate, method = "spearman")




####3rd vaccination at the end of november
data9 <- merge(data, vacc3_nov)
data9$ex_diff_q0.5_2021_month <- data9$ex_diff_q0.5_2021 * 12


data9 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Boosterd 60+ at the end of November", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/vacc3_november.png", plot = last_plot(), dpi = 100)


##regression modell
model_60_nov_20_21_3rd <- lm(ex_diff_q0.5_2021_month ~ rate, data=data9)
summary(model_60_nov_20_21_3rd)


###correlations
cor.test(data9$ex_diff_q0.5_2021_month, data9$rate, method = "pearson")
cor.test(data9$ex_diff_q0.5_2021_month, data9$rate, method = "spearman")






























































####time until 30 % got vaccinated
data5 <- merge(data, min30)
data5$ex_diff_q0.5_2021_month <- data5$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data5, aes(Date, ex_diff_q0.5_2021_month, color = country_name)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data5 %>% 
  ggplot() +
  geom_point(mapping = aes(x = Date, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Time where 30% of Population 60+ were Vaccinated", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Timing of Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_30perc_60.png", plot = last_plot(), dpi = 100)

##regression modell
model <- lm(ex_diff_q0.5_2021_month ~ Date, data=data5)
summary(model)


###correlations
#cor.test(data5$ex_diff_q0.5_2021_month, data5$Date, method = "pearson")



####time until 50% vaccinated
data6 <- merge(data, min50)
data6$ex_diff_q0.5_2021_month <- data6$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data6, aes(Date, ex_diff_q0.5_2021_month, color = country_name)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data6 %>% 
  ggplot() +
  geom_point(mapping = aes(x = Date, y = ex_diff_q0.5_2021_month, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Time where 50% of Population 60+ were Vaccinated", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Timing of Vaccination Uptake 60+")
ggsave("U:/gits/ex2021/tmp/vacc_50perc_60.png", plot = last_plot(), dpi = 100)
##regression modell
model <- lm(ex_diff_q0.5_2021_month ~ Date, data=data6)
summary(model)


###correlations
cor.test(data6$ex_diff_q0.5_2021_month, data6$rate, method = "pearson")



####vacc until august for countries with lost

data8 <- merge(data, vacc2_august_lost)
data8$diff_all <- (data8$ex_diff_q0.5_2020 + data8$ex_diff_q0.5_2021)*12
data8$ex_diff_q0.5_2021_month <- data8$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data8, aes(rate, ex_diff_q0.5_2021_month, color = country_name)) + 
  geom_point() + 
  geom_smooth(method = "lm")

data8 %>% 
  ggplot() +
  geom_point(mapping = aes(x = rate, y = diff_all, color = Country, shape = Country)) +
  scale_shape_manual(values=1:27) +
  labs(x="Percentage Vaccinated at the end of August", y="Change in e0 from 2020 to 2021 in months", 
       title = "Change in e0 from 2020 to 2021 and Vaccination Uptake")
ggsave("U:/gits/ex2021/tmp/total_vacc_september_lost.png", plot = last_plot(), dpi = 100)

##regression modell
model <- lm(ex_diff_q0.5_2021_month ~ rate, data=data8)
summary(model)


###correlations
cor.test(data8$ex_diff_q0.5_2021_month, data8$rate, method = "pearson")


