library(dplyr)
rates <- read_rds("U:/gits/ex2021/dat/coverage/rates.rds") %>% 
  filter(Country != "Lithuania")
totals.rates <- read_rds("U:/Jonas/totals_rates.rds")%>% 
  filter(Country != "Lithuania")
###fill in missing dates
vacc2 <- rates %>% 
  filter(Measure == "Vaccination2") %>% 
  group_by(Country) %>% 
  mutate(Date = as.Date(Date)) %>%
  tidyr::complete(Date = seq.Date(min(Date), max(Date), by="day")) %>% 
  fill(`rate`)

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

min30 <- vacc2 %>% 
  filter(rate >=0.3) %>% 
  group_by(Country) %>% 
mutate(ticker = row_number()) %>% 
  filter(ticker == "1",
         Country != "Finland",
         Country != "Northern Ireland")

min50 <- vacc2 %>% 
  filter(rate >=0.5) %>% 
  group_by(Country) %>% 
  mutate(ticker = row_number()) %>% 
  filter(ticker == "1",
         Country != "Finland",
         Country != "Northern Ireland")


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


data <- as_tibble(fig$e0diff$data)
names(data)[25] <- "Country"
data2 <- merge(data, vacc2_august)
data2$diff_all <- (data2$ex_diff_q0.5_2020 + data2$ex_diff_q0.5_2021)*12
data2$ex_diff_q0.5_2021_month <- data2$ex_diff_q0.5_2021 * 12

##plotting with linear fit line
ggplot(data2, aes(rate, diff_all, color = Country)) + 
  geom_point() + 
  geom_smooth(method = "lm")

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
##regression modell
model1 <- lm(diff_all ~ rate, data=data2)
summary(model1)

model2 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data2)
summary(model1)
###correlations
cor.test(data2$ex_diff_q0.5_2021_month, data2$rate, method = "pearson")



####vaccination at the end of september
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
model2 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data3)
summary(model2)


###correlations
cor.test(data3$ex_diff_q0.5_2021_month, data3$rate, method = "pearson")



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
model3 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data4)
summary(model3)


###correlations
cor.test(data4$ex_diff_q0.5_2021_month, data4$rate, method = "pearson")


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
model3 <- lm(ex_diff_q0.5_2021_month ~ rate, data=data7)
summary(model3)


###correlations
cor.test(data7$ex_diff_q0.5_2021_month, data7$rate, method = "pearson")


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


