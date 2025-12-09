################################################################
################################################################
#🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃#
#🎃                                                          🎃#
#🎃                                                          🎃#
#🎃     This is my R Script Submission for the Pumpkins      🎃#
#🎃                      Dataset.                            🎃#
#🎃          Author: Shahwar Nadeem, 20430461                🎃#
#🎃             Module Code: LIFE4138_2526                   🎃#
#🎃                                                          🎃#
#🎃                                                          🎃#
#🎃                                                          🎃#
#🎃                                                          🎃#
#🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃🎃#
################################################################
################################################################
library(ggplot2)
library(tidyverse)
library(dplyr)
library(magrittr)
library(readr) 
pumpkins_12 <- read_csv("pumpkins_12.csv")
View(pumpkins_12)

pumpkins_12 %>%
  arrange(desc(weight_lbs))
#Code to return maximum Pumpkin weight

pumpkins_12$Weight_kg <- pumpkins_12$weight_lbs * 0.45359237
#Create a new column titled "Weight_kg" which will convert column "weight_lbs" into kg

pumpkins_12 %>%
  select(id, variety, weight_lbs, Weight_kg)
#Display the id, variety, weight_lbs, weight_kg columns for pumpkins

pumpkins_12$weight_class <- ifelse (pumpkins_12$weight_lbs >= 1000, "Heavy",
                                    ifelse (pumpkins_12$weight_lbs >= 200, "Medium", "Light"))
#ifelse statements to return strings "Heavy", "Medium" or "Light" depending on weight_lbs

ggplot(aes(x = est_weight, y = weight_lbs, colour = weight_class), data = pumpkins_12)+
  geom_point(na.rm = TRUE)+
  labs(title = "Estimated Weight vs Actual Weight of Pumpkins",
       x = "Estimated Weight (lbs)",
       y = "Weight (lbs)")
#Create a Scatter Plot displaying Estimated vs Actual Pumpkin Weight, colour by class

pumpkins_filtered <- pumpkins_12 %>%
  filter( country == "Italy" | country == "Germany" | country == "France") %>%
  select(id, variety, country, Weight_kg, weight_class)
write.csv(pumpkins_filtered, "pumpkins_filtered.csv")
View(pumpkins_filtered)
#Filter dataset for Pumpkins from France, Italy and Germany only, Save filtered dataset as pumpkins_filtered.csv file

pumpkins_12 %>%
  group_by(country, variety) %>%
  summarise(mean_weight = mean(Weight_kg, na.rm = TRUE)) %>%
              select(variety, country, mean_weight)
#Mean weight for each variety of pumpkin from each country in original dataset

pumpkins_filtered %>%
  group_by(country, variety) %>%
  summarise(mean_weight = mean(Weight_kg, na.rm = TRUE)) %>%
  select(variety, country, mean_weight) %>%
  arrange(mean_weight)
#Mean weight for each variety of pumpkin from Italy, Germany, France in filtered dataset, beginning with the smallest mean weight

pumpkins_filtered %>%
  group_by(country) %>%
  summarise(mean_weight = mean(Weight_kg, na.rm = TRUE)) %>%
  select(country, mean_weight)
#Mean weight for pumpkins from Italy, Germany, France in filtered dataset

ggplot(aes(country, Weight_kg), data = pumpkins_filtered)+
         geom_boxplot()+
  labs(title = "Distribution of Pumpkin Weight between Countries",
          x = "Country",
          y = "Weight (kg)")
#Boxplot of the pumpkin weight for the countries in filtered dataset

ggplot(aes(country, Weight_kg, colour = country), data = pumpkins_filtered)+
  geom_boxplot()+
  facet_wrap(~ variety)+
  labs(title = "Distribution of Pumpkin Weight between Varieties in Countries",
       x = "Country",
       y = "Weight (kg)")
#Boxplot of the pumpkin weight for countries in filtered dataset for each variety of pumpkin