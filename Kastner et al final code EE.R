### Kastner et al. manuscript analysis and graphing code
### Submitted to Ecology and Evolution 13MAR2024

### load libraries  --------

library(readr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(lubridate)
library(aod)
library(boot)
library(ggResidpanel)

### load data --------

mortality <- read_excel("data/Kastner et al data for submission EE.xlsx", sheet = 1, col_types = NULL, col_names =T, na = "NA")
banding <- read_excel("data/Kastner et al data for submission EE.xlsx", sheet = 2, col_types = NULL, col_names =T, na = "NA")
aafb_snakes <- read_excel("data/Kastner et al data for submission EE.xlsx", sheet = 3, col_types = NULL, col_names =T, na = "NA")
ebl_snakes <- read_excel("data/Kastner et al data for submission EE.xlsx", sheet = 4, col_types = NULL, col_names =T, na = "NA")
rpm_snakes <- read_excel("data/Kastner et al data for submission EE.xlsx", sheet = 5, col_types = NULL, col_names =T, na = "NA")

### 1. Question 1: frequency of unsuccessful ingestion  -------- 

# 1a. check column names and structure 

colnames(mortality)

str(mortality)

# 1b. remove NAs, and also any unbanded individuals 

mortality <- mortality %>%
  filter(complete.cases(cod,band_num,tx))

# 1c. filter out 2017-2018 birds 

D <- as.Date("2019-11-24")

mortality$band_date <- ymd(mortality$band_date)

mortality <- mortality %>%
  filter(mortality$band_date > D)

mortality <- mortality %>%
  mutate_at(vars(cod, name, sex), factor) 

summary(mortality$cod)

# 1d. make a new column for x axis

# goal is to have bts, cat, other and unknown as 4 categories on x axis, and stack 
# slimed (not successfully ingested) and consumed

mortality <- mortality %>%
  mutate(pred = case_when(
    cod == "slimed by bts" ~ "bts",
    cod == "bts predation" ~ "bts",
    cod == "cat predation" ~ "cat",
    cod == "other" ~ "other",
    cod == "unknown" ~ "unknown"))

pred_order <- factor(mortality$pred, levels = c('bts', 'cat', 'other', 'unknown'))

# 1e. graph ----

# png("figures/figure 2.png", width = 5000, height = 2800,
#     units = "px", res = 450, bg = "white")

ggplot(mortality, aes(x = factor(pred_order)))+
  geom_bar_pattern(stat = "count", color = "black", pattern = "none", width = 0.55, aes(fill = cod)) + # remove pattern = "none" & add pattern_type = cod inside the aes?
  # scale_pattern_type_discrete(values = c("slimed by bts" = "stripe",
  #   "bts predation" = none,
  #   "cat predation" = none,
  #   "other" = none,
  #   "unknown" = none)) +
  theme_classic() +
  scale_fill_manual(values = c("gray87", "gray87", "gray87", "white", "gray87")) +
  theme(legend.position = "none") +
  xlab("Cause of mortality") + 
  ylab("Number of records") +
  scale_x_discrete(labels = c("bts" = "BTS predation", 
                              "cat" = "Cat predation", 
                              "other" = "Other causes",
                              "unknown" = "Unknown causes")) +
  annotate("text", x = 1, y = 65, label = "Slimed", size = 6) +
  annotate("text", x = 1, y = 138, label = "Consumed", size = 6)+
  theme(text = element_text(size = 20)) 

# dev.off()

### 2. Question 2: size of prey vs ingestion success    --------

# 2a. merge banding and mortality databases

banding <- banding %>%
  filter(!is.na(tx_num)) # remove entries with no transmitter number

band_mort <- mortality %>%
  left_join(banding, by = c("band_num"))

# 2b. only consider cases where BTS were cause of death

band_mort <- band_mort %>%
  filter(cod == "bts predation"|cod == "slimed by bts") %>%
  filter(complete.cases(cod))

# 2c. format dates and remove entries from 2017-2018

names(band_mort)

D <- as.Date("2019-11-24")

band_mort$band_date <- ymd(band_mort$band_date)

band_mort <- band_mort %>%
  filter(band_mort$band_date > D)

# 2d. make a new column with binary cod

band_mort <- band_mort %>%
  mutate(binary_cod = if_else(cod == "slimed by bts",1,0))

# names(band_mort)

# 2e. run logistic regression 

model <- glm(formula = binary_cod ~ mass, family = binomial (link = 'logit'), data = band_mort)

resid_panel(model)
plot(model$residuals, pch = 16, col = "red") #check model fit

#results

summary(model) # z value = -0.515, p = 0.607 - Wald test

# 2f. graph

# ggplot(band_mort, aes(x= mass, y= binary_cod)) + geom_point() +
#   stat_smooth(method="glm", color="green", se=FALSE,
#               method.args = list(family=binomial))+
#   theme_classic()+
#   scale_y_continuous(breaks = seq(0,1,1), labels=c("0" = "Consumed", 
#                                                    "1" = "Slimed"))+
#   
#   xlab("Fledgling mass (g)")+
#   ylab("Fate")+
#   theme(text = element_text(size = 20),
#         axis.text.y = element_text(angle=45,
#                                    hjust = 0.5))


### 3. Question 3: size of pred vs ingestion success    --------

# 3a. take subset that eat endothermic prey (larger than 800 svl)

aafbsnakes_endo <- aafb_snakes %>%
  subset(svl>900)

svl_aafb <- aafbsnakes_endo$svl # pull out a vector of svl's (snout vent lengths) only

# 3b. resample 52 snakes out of the overall sample 5000 times

bootmeans <- tibble(num = 1:5000) %>% 
  group_by(num) %>%
  mutate(means = mean(sample(svl_aafb, size = 52, replace = TRUE)))

quant_snakes <- quantile(bootmeans$means, c(0.025,0.975)) # quantile means

print(quant_snakes)


# 3c. graph

bootmeansmean <- mean(bootmeans$means) #mean of bootstrap means
bootmeansmean

# merge ebl snakes (successful ingestion cases) with visual survey sample (aafb snakes)

all_snakes <- ebl_snakes %>%
  full_join(aafb_snakes, (by = c("pit", "clip", "svl", "tl", "mass", "sex", "id")))

all_snakes<-all_snakes%>%
  distinct(pit, .keep_all = T) # remove duplicates 

all_snakes<-all_snakes%>%
  filter(complete.cases(svl)) # remove incomplete rows for svl

all_snakes <- all_snakes %>%
  mutate(id2 = case_when(
    id == "ebl" ~ "ebl",
    id == "aafb" & svl < 900 ~ "aafb < 900",
    id == "aafb" & svl >= 900 ~ "aafb > = 900" 
  )) # recode based on values in svl column

# 3d. graph

# png("figures/figure 3.png", width = 5000, height = 2800,
#     units = "px", res = 300, bg = "white")

ggplot(all_snakes, aes(x = svl, fill = id2)) +                       # Draw overlaying histogram
  geom_histogram(aes(y = stat(count/sum(count))), position = "identity", alpha = 0.4, bins = sqrt(nrow(all_snakes)))+
  theme_classic()+
  xlab("SVL (mm)")+
  ylab("Proportion")+
  geom_segment(x = 1114, xend = 1114, y = 0, yend = 0.17, alpha = 0.5, color = "black", linetype = "dashed") + #mean size bootstrapped
  # annotate("rect", xmin = 1068, xmax = 1167, ymin = 0, ymax = 0.1149, alpha = 0.1, fill = "gray4") + # 2.5 and 97.5% quantile means
  # annotate("rect", xmin = 1068, xmax = 1167, ymin = 0.1149, ymax = 0.17, alpha = 0.4, fill = "gray4") +
  geom_segment(x = 1206, xend = 1206, y = 0, yend = 0.17, alpha = 0.5, color = "black", size = 0.8) + #mean size confirmed ingestion
  scale_fill_manual(values = c("darkgray", "skyblue", "navy blue"), name = "Sample", labels = c("Visual search (< 800 mm SVL)", "Visual search (≥ 800 mm SVL)", "Consumed Såli"))+
  theme(text = element_text(size = 24)) 

# dev.off()

### 4. Question 4: RPM for successful attempts          --------

# 4a. Rename columns

rpm_snakes <- rpm_snakes %>%
  rename(prey_mass = 'mass.y',
         pred_mass = 'mass.x')

# 4b. Calculate RPM

rpm_snakes <- rpm_snakes %>%
  mutate(rpm_ratio = prey_mass / pred_mass,
         rpm_perc = (prey_mass/pred_mass) * 100)

# 4b. Remove last 2 rows because they are two birds that were consumed by the same snake

rpm_snakes <- rpm_snakes %>%
  slice_head(n = 54)

# 4c. Linear model

lm_rpm <- lm(prey_mass ~ pred_mass, data = rpm_snakes)

summary(lm_rpm) # t = -0.2, p = 0.604

summary(lm_rpm)$r.squared

resid_panel(lm_rpm)

# 4d. graph

# png("figures/figure 4.png", width = 5000, height = 2800,
#     units = "px", res = 400, bg = "white")

ggplot(rpm_snakes, aes(x = pred_mass, y = prey_mass)) +
  geom_point(position=position_jitter(h=1,w=1), size = 2) +
  # geom_smooth(method = lm, colour = "grey27", size = 0.5) +
  theme_classic() +
  xlab("Snake mass (g)") +
  ylab("Fledgling mass (g)")+
  scale_y_continuous(limit = c(55, 90))+
  theme(text = element_text(size = 20)) +
  geom_abline(slope = 0.2, size = 0.5, colour = "gray50") +
  annotate("text", x = 470, y = 90, label = "20%", size = 5, colour = "grey30") +
  geom_abline(slope = 0.4, size = 0.5, colour = "gray50") +
  annotate("text", x = 245, y = 90, label = "40%", size = 5, colour = "grey30") +
  geom_abline(slope = 0.6, size = 0.5, colour = "gray50") +
  annotate("text", x = 170, y = 90, label = "60%", size = 5, colour = "grey30") +
  geom_abline(slope = 0.8, size = 0.5, colour = "gray50") +
  annotate("text", x = 130, y = 90, label = "80%", size = 5, colour = "grey30") 

# dev.off()