#load packages####

library(Hmisc) # for correlation matrix
library(visreg)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(emmeans)
library(lubridate)
library(lmerTest)
library(factoextra)
library(gamm4)


#load data####

maze_df <- read.csv('summary_maze.csv', fileEncoding="UTF-8-BOM")
control_df <- read.csv('summary_control.csv', fileEncoding="UTF-8-BOM")
solving.time_df <- read.csv('solving_time.csv', fileEncoding="UTF-8-BOM")
combined_df <- read.csv('combined_df.csv', fileEncoding="UTF-8-BOM")

solving.time_df$level <- as.factor(solving.time_df$level)

maze_df <- maze_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 9),
         cos_cycle = cos(2 * pi * cycle / 9),)
control_df <- control_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 9),
         cos_cycle = cos(2 * pi * cycle / 9),)
control_df2 <- control_df2 %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 9),
         cos_cycle = cos(2 * pi * cycle / 9),)
solving.time_df <- solving.time_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 9),
         cos_cycle = cos(2 * pi * cycle / 9),)
combined_df <- combined_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 9),
         cos_cycle = cos(2 * pi * cycle / 9),)


#PCA####

##PCA:
###perform PCA:

behavior_matrix <- combined_df %>% select(14:19)
behavior_matrix <- scale(behavior_matrix)
PCA_results <- prcomp(behavior_matrix, scale = TRUE)

PCA_loadings <- PCA_results$rotation
PCA_eigenvalues <- get_eigenvalue(PCA_results)

###plot variables' contribution to the first three PCs:

fviz_contrib(PCA_results, choice = "var", axes = 1) +
  labs(title = "Contributors to PC1",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) + 
  theme( axis.text.x  = element_text( angle = 45, hjust = 1)) 

fviz_contrib(PCA_results, choice = "var", axes = 2, top = 10) +
  labs(title = "Contributors to PC2",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) + 
  theme( axis.text.x  = element_text( angle = 45, hjust = 1)) 

fviz_contrib(PCA_results, choice = "var", axes = 3, top = 10) +
  labs(title = "Contributors to PC3",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) + 
  theme( axis.text.x  = element_text( angle = 45, hjust = 1)) 


###PCA details and stats:
####creating df:
PCA_loadings <- PCA_results$rotation
scores <- as.data.frame(PCA_results$x[,1:3]) %>%
  mutate(sin_cycle = combined_df$sin_cycle, 
         cos_cycle = combined_df$cos_cycle, 
         cycle = combined_df$cycle, 
         species = combined_df$species,
         treatment = combined_df$treatment)
PCA_eigenvalues <- get_eigenvalue(PCA_results)

scores.plot.maze <- scores[which(scores$treatment == 'maze'),] %>%
  bind_rows(scores[which(scores$treatment == 'maze'),] %>%
              filter(cycle == "1") %>%         
              mutate(cycle = 9))

scores.plot.control <- scores[which(scores$treatment == 'control'),] %>%
  bind_rows(scores[which(scores$treatment == 'control'),] %>%
              filter(cycle == "1") %>%         
              mutate(cycle = 9))

write.csv(PCA_loadings,'PCA_loadings.csv')

write.csv(PCA_eigenvalues,'PCA_eigenvalues.csv')
#PC1####
##model PC1
model <- lm(PC1 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC1 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC1 ~ cos_cycle*species*treatment , 
                data = scores)


gam_model <- gamm4(
  PC1 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores
)


pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor(pred_gam, pred_m2) #best model
cor(pred_gam, pred_m3)
cor(pred_gam, pred_m4)

summary(cos_model)

cos_model_2 <- lm(PC1 ~ cos_cycle+treatment+species, 
                  data = scores)

summary(cos_model_2)

model_PC1.maze <- summary(cos_model_2)
write.csv(model_PC1.maze$coefficients, 'model_PC1.csv')


##plot PC1
ggplot(scores, 
       aes(x = cycle, y = PC1)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ cos(2 * pi * x/9),
              se = TRUE)   + 
  labs(x = "Treatment",y = "PC1",
       title = "PC1: foraging related activities") +
  theme_classic(base_size = 25) +
  scale_x_continuous(breaks = c(5),
                     labels = c("Full Moon"))  + 
  theme(legend.position = "none")+ 
  facet_wrap(~treatment)


#PC2####
##model PC2
model <- lm(PC2 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC2 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC2 ~ cos_cycle*species*treatment , 
                data = scores)


gam_model <- gamm4(
  PC2 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores
)


pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor(pred_gam, pred_m2) #best model
cor(pred_gam, pred_m3)
cor(pred_gam, pred_m4)

summary(cos_model)

cos_model_2 <- lm(PC2 ~ cos_cycle + species + treatment , 
                  data = scores)

summary(cos_model_2)

model_PC2.maze <- summary(cos_model_2)
write.csv(model_PC2.maze$coefficients, 'model_PC2.csv')


##plot PC2 - maze
ggplot(scores, 
       aes(x = cycle, y = PC2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ cos(2 * pi * x/9),
              se = TRUE)   + 
  labs(x = "Moon cycle",y = "PC2",
       title = "PC2: vigilance & digging") +
  theme_classic(base_size = 25)  +
  scale_x_continuous(breaks = c(5),
                     labels = c("Full Moon"))  + 
  theme(legend.position = "none")+ 
  facet_wrap(~treatment)



#PC3####
##model PC3
model <- lm(PC3 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC3 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC3 ~ cos_cycle*species*treatment , 
                data = scores)

gam_model <- gamm4(
  PC3 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores
)


pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor(pred_gam, pred_m2) #best model
cor(pred_gam, pred_m3)
cor(pred_gam, pred_m4)

summary(sin_model)

sin_model_2 <- lm(PC3 ~ sin_cycle + species + treatment , 
                  data = scores)

summary(sin_model_2)

model_PC3.maze <- summary(sin_model_2)
write.csv(model_PC3.maze$coefficients, 'model_PC3.csv')


##plot PC3 - maze
ggplot(scores, 
       aes(x = cycle, y = PC3)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ sin(2 * pi * x/9),
              se = TRUE)   + 
  labs(x = "Moon cycle",y = "PC3",
       title = "PC3: competition") +
  theme_classic(base_size = 25) +
  scale_x_continuous(breaks = c(5),
                     labels = c("Full Moon"))  + 
  theme(legend.position = "none")+ 
  facet_wrap(~treatment)


#correlations between variables####
cols <- c("second_tray", "combined.seeds", "explore_score")
cor_matrix <- as.matrix(maze_df[which(maze_df$trial_length == 60), cols])

res <- rcorr(cor_matrix, 
             type = "pearson")

res$r   # correlations
res$P   # p-values

#overall seeds taken####

seeds_df <- combined_df %>% 
  filter(trial_length == 60) %>%
  mutate(overall_seeds = seeds_maze + seeds_control) %>%
  select(treatment, species, cycle, sin_cycle, cos_cycle, overall_seeds, 
         seeds_maze, seeds_control) 

#model maze
model <- lm(overall_seeds ~ sin_cycle * species * treatment +
              cos_cycle * species * treatment, 
            data = seeds_df)
sin_model <- lm(overall_seeds ~ sin_cycle * species * treatment, 
                data = seeds_df)
cos_model <- lm(overall_seeds ~ cos_cycle * species * treatment, 
                data = seeds_df)

gam_model <- gamm4(
  overall_seeds ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = seeds_df
)


pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor.test(pred_gam, pred_m2) #best model
cor.test(pred_gam, pred_m3)
cor.test(pred_gam, pred_m4)


summary(model)

model2 <- lm(overall_seeds ~ sin_cycle + species + treatment +
               cos_cycle, 
             data = seeds_df)

summary(model2)

model_seeds_maze <- summary(model2)
write.csv(model_seeds_maze$coefficients, 'model_seeds.csv')


#plot

seeds_df.plot <- seeds_df %>%
  bind_rows(seeds_df %>%
              filter(cycle == "1") %>%         
              mutate(cycle = 9))

ggplot(seeds_df, aes(x = cycle, y = seeds_control)) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ sin(2 * pi * x/9) + cos(2 * pi * x/9),
    se = TRUE)  +
  labs(x = "Moon phase", y = "Seeds taken (g)", 
       title = "Overall seeds taken from trays") +
  theme_classic(base_size = 25) +
  scale_x_continuous(
    breaks = c(5),
    labels = c("Full Moon")
  )  + theme(legend.position = "none") + 
  facet_wrap(~treatment)


#time to find the second tray####


##model selection - control
model <- lm(second_tray ~ sin_cycle * species * treatment +
              cos_cycle * species * treatment, 
            data = combined_df)
sin_model <- lm(second_tray ~ sin_cycle * species * treatment, 
                data = combined_df)
cos_model <- lm(second_tray ~ cos_cycle * species * treatment, 
                data = combined_df)

gam_model <- gamm4(
  second_tray ~
    s(cycle, by = interaction(treatment, species), k = 4) +
    treatment * species,
  data = combined_df
)


pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor.test(pred_gam, pred_m2) #best model
cor.test(pred_gam, pred_m3)
cor.test(pred_gam, pred_m4)

log_model <- lm(log(second_tray) ~ sin_cycle * species * treatment +
                  cos_cycle * species * treatment, 
                data = combined_df)

AIC(model, log_model)

summary(log_model)

log_model2 <- lm(log(second_tray) ~ sin_cycle * treatment +
                   cos_cycle * treatment + species, 
                 data = combined_df)

summary(log_model2)

plot(log_model2)

model_second <- summary(log_model2)
write.csv(model_second$coefficients, 'model_second.csv')


##plot
ggplot(combined_df, aes(x = cycle, y = log(second_tray))) +
  geom_point(alpha = 0.6) + 
  ylim(-2,5) +
  geom_smooth(
    method = "lm",
    formula = y ~ sin(2 * pi * x/9) + cos(2 * pi * x/9),
    se = TRUE)  +
  labs(x = "Moon phase", y = "Time to discover second tray (log[min])", 
       title = "Time to discover second tray") +
  theme_classic(base_size = 25) +
  scale_x_continuous(
    breaks = c( 5),
    labels = c( "Full Moon")
  ) + theme(legend.position = "none") + 
  facet_wrap(~treatment)
