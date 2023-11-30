## load packages ##

library(readxl)
library(dplyr)
library(agricolae)
library(FSA)


## Code for Kruskal-Wallis and Dunn's tests applied to Figure 2.

library(readxl)


## load mating system data ##
matingsystem <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/Callen_matingsystem_foranalysis_4July2023.xlsx")


## eliminate rows that do not contain male and female interactions ##
matingsystem_foranalysis <- matingsystem[!is.na(matingsystem$Number_of_males_and_females),]


## Kruskal test and Dunn's tests for Fast chase rate across species (Fig. 2A)

kruskal_fastchase <- kruskal.test(Fast_chases_per_minute ~ Species, data = matingsystem_foranalysis)

kruskal_fastchase

dunn_fastchase <- dunnTest(Fast_chases_per_minute ~ Species,
                            data=mating_system_foranalysis,
                            method="bonferroni")

dunn_fastchase

## Kruskal test and Dunn's tests for Time in frontal display per minute across species ##


kruskal_display <- kruskal.test(Time_frontal_display_minute ~ Species, data = matingsystem_foranalysis)

kruskal_display

dunn_display <- dunnTest(Time_frontal_display_minute ~ Species,
                           data=mating_system_foranalysis,
                           method="bonferroni")

dunn_display




## Code for analysis reported in Figure 3 ##



## load dataframes ##

Sociability <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/Sociability_MasterSheet_27Oct.xlsx")
## this loads dataframe containing sociability data across species, used in Figure 3A and 3C ##

## Kruskal test and Dunn's tests for Sociability (Fig. 3C) ##

kruskal_sociability <- kruskal.test(TimeinInteraction ~ Species, data = Sociability)

kruskal_sociability

dunn_sociability <- dunnTest(TimeinInteraction ~ Species,
                           data=Sociability,
                           method="bonferroni")

dunn_sociability

## Kruskal test and Dunn's tests for activity in the shoaling/sociability assay ##


kruskal_activity <- kruskal.test(Sociability_distance_moved ~ Species, data = Sociability)

kruskal_activity

dunn_activity <- dunnTest(Sociability_distance_moved ~ Species,
                         data=Sociability,
                         method="bonferroni")

dunn_activity

## Load scototaxis dataframe. These data are used in Figure 3B and 3D ##

Scototaxis <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/Scototaxis_30Oct2023.xlsx")

## Kruskal test and Dunn's tests for stress movement in the scototaxis assay, Fig. 3B ##


kruskal_stressmovement <- kruskal.test(Scototaxis_distance_moved ~ Species, data = Scototaxis)

kruskal_stressmovement

dunn_stressmovement <- dunnTest(Sociability_distance_moved ~ Species,
                          data=Sociability,
                          method="bonferroni")

## Kruskal test and Dunn's tests for boldness in the scototaxis assay, Fig. 3D ##


kruskal_boldness <- kruskal.test(TimeinWhite ~ Species, data = Scototaxis)

kruskal_boldness

dunn_boldness <- dunnTest(TimeinWhite ~ Species,
                                data=Scototaxis,
                                method="bonferroni")


## load route learning dataframe containing route learning performance across trials (Fig. 3F) ##
RouteLearning <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/RL_full_31Oct2023.xlsx")



## run linear mixed-effects model telling us the effect of species and trial on number of wrong door entries per trial ##
error_rate_model <- lmer(Wrong_door_entries ~ Species + trial + (1 | fishID),
              data = RouteLearning)
summary(error_rate_model)

## Determine slope of wrong door entries between trials 1 and 3 for each species and whether or not significantly different from 0 ##


library(dplyr)

## filter out trials 4 and 5 ##

RouteLearning1_3 <- filter(RouteLearning, trial == "1" | trial == "2" | trial == "3")

## build regression models ##

regression_models = RouteLearning %>% group_by(Species) %>% do(model = lm(wrong_door_entries ~ trial, data = .))

## extract slope and associated p-value ##

summary(regression_models)

## load Detour data sheet for Figs. 3E and 3G

Detour <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/Detour_31Oct2023.xlsx")

## Kruskal test and Dunn's tests for Detour social motivation (Fig. 3E) ##

kruskal_socialmotivation <- kruskal.test(Detour_timetoglass ~ Species, data = Detour)

kruskal_socialmotivation

dunn_socialmotivation <- dunnTest(Detour_timetoglass ~ Species,
                            data=Detour,
                            method="bonferroni")

dunn_socialmotivation

## Kruskal test and Dunn's tests for Detour solve time (Fig. 3G) ##

kruskal_detour <- kruskal.test(Detour_solve_time ~ Species, data = Detour)

kruskal_detour

dunn_detour <- dunnTest(Detour_timetoglass ~ Species,
                                      data=Detour,
                                      method="bonferroni")

dunn_detour


## this loads dataframe containing the total number of wrong door entries summed across trials ##

RL_error_rate <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/RL_2021_for_multivariate.xlsx")



kruskal_errors <- kruskal.test(Total_wrong_door_entries ~ Species, data = RL_error_rate)

kruskal_errors

dunn_errors <- dunnTest(Total_wrong_door_entries ~ Species,
                            data=RL_error_rate,
                            method="bonferroni")

dunn_errors




## Below is the code used in the Linear Discriminant Analysis, the results of which are presented in Figure 4 ##


## load complete dataset (excluding L. perugiae) ##
FullAssayData <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/LDA_combined_sheet_14Sep.xlsx")

FullAssayData %>%
  select(Species, Detour_solve_time, TimeinWhite, Scototaxis_distance_moved, TimeinInteraction, Sociability_distance_moved, error_rate_slope_1to3) %>% # make sure we select only the columns we care about
  mutate(Species = factor(Species, levels = c("X. nigrensis", "H. formosa", "G. vittata"))) # make sure we're consistent about ordering (complex, then simple)


## Construct LD1 and LD2 calculate contributions of each one to between-species variation in behavior and cognition ##

r <- lda(formula = Species ~ .,
         data = FullAssayData,
         prior = c(1,1,1)/3)

prop = r$svd^2/sum(r$svd^2)

prop

## add LD axes to dataframe ##

lda_fit <- FullAssayData %>%
  mutate_at(vars(error_rate_slope_1to3:TimeinInteraction), list(~ scale(., center = TRUE))) %>% # center and scale each variable
  lda(Species ~ Sociability_distance_moved + error_rate_slope_1to3 + Detour_solve_time + TimeinWhite + Scototaxis_distance_moved + TimeinInteraction + Sociability_distance_moved, data = .)

## add the loadings back to the dataframe ##
transform <- predict(lda_fit)$x
df_LDA <- data_for_lda %>%
  mutate(
    lda1 = transform[, 1],
    lda2 = transform[, 2]
  )

## compare species in LD1 ##
anova <- aov(lda1 ~ Species, data = df_LDA)

summary(anova)

TukeyHSD(anova)

## compare species in LD2 ##

anova2 <- aov(lda2 ~ Species, data = df_LDA)

summary(anova2)

TukeyHSD(anova2)

## Code used to generate SEM results reported in Figure 5.##
## Note: I demonstrate this with one of the two species with scaled data, and I show a single example of an SEM. I ran 57 models for each species, and showing the code for those would require too many lines. ##

## load lavaan package ##
library(broom)
library(lavaan)

## load dataframe with just G. vittata behavior and cognition ##

vittata <- read_excel("C://Users/scott/OneDrive - The University of Texas at Austin/DataSheets/vittata_Callen_fourassays_14Sep.xlsx")

cols <- c("TimeinWhite", "TimeinInteraction", "Detour_solve_time", "error_rate_slope_1to3", "Sociability_distance_moved", "Scototaxis_distance_moved")

vittata_2 <- vittata[cols]


j <- sapply(vittata_2, is.numeric)
## scale data in data frame. Note: this was only done for H. formosa and G. vittata, but not X. nigrensis ##
vittata_2[j] <- scale(vittata_2[j])

## model specification. This is a model with all 6-variables ##

model1 <- 'Syndrome =~ Scototaxis_distance_moved + Detour_solve_time + error_rate_slope_1to3 + Sociability_distance_moved'

## run the SEM ##
sem.fit.full1 <- sem(model1, data = nigrensis_2, check.gradient = FALSE, std.lv = TRUE, missing = "pairwise")
summary(sem.fit.full1, standardized = TRUE, fit.measures = TRUE)

## extract fit measures (GFI, RMSEA, Chi-squared p-value) to determine whether model meets fit criteria for AIC comparison ##
fitMeasures(sem.fit.full1)




