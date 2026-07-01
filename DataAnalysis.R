#load packages####

library(Hmisc) # for correlation matrix
library(ggplot2) #for plotting
library(tidyverse) #for data wrangling  
library(dplyr) #for data wrangling  
library(emmeans) #for post-hoc analyses 
library(lubridate) #for time and date data
library(lmerTest) #for linear models
library(factoextra) #for PCA
library(gamm4) #for model selection
library(ggeffects) #for plotting

#load data####

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load data
maze_df <- read.csv('summary_maze.csv', fileEncoding="UTF-8-BOM")
control_df <- read.csv('summary_control.csv', fileEncoding="UTF-8-BOM")
combined_df <- read.csv('combined_df.csv', fileEncoding="UTF-8-BOM")

#transform lunar cycle into sine and cosine waves
maze_df <- maze_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 8),
         cos_cycle = cos(2 * pi * cycle / 8),)
control_df <- control_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 8),
         cos_cycle = cos(2 * pi * cycle / 8),)
combined_df <- combined_df %>%
  mutate(sin_cycle = sin(2 * pi * cycle / 8),
         cos_cycle = cos(2 * pi * cycle / 8),)


#PCA####

##perform PCA####

#create a behavior matrix
behavior_matrix <- combined_df %>% select(14:19)

#scale
behavior_matrix <- scale(behavior_matrix)

#perform PCA
PCA_results <- prcomp(behavior_matrix, scale = TRUE)

#save results
PCA_loadings <- PCA_results$rotation

PCA_eigenvalues <- get_eigenvalue(PCA_results)

scores <- as.data.frame(PCA_results$x[,1:3]) %>%
  mutate(sin_cycle = combined_df$sin_cycle, 
         cos_cycle = combined_df$cos_cycle, 
         cycle = combined_df$cycle, 
         species = combined_df$species,
         treatment = combined_df$treatment)
scores$treatment <- as.factor(scores$treatment)


##plot variables' contribution to the first three PCs####

#order behaviors in a fixed order
my_order <- c(
  "f_maze",
  "f_control",
  "open",
  "competition",
  "dig",
  "antipredator")

my_labels <- c(
  f_maze = "forage tray 1",
  f_control = "forage tray 2",
  open = "open",
  competition = "aggresive",
  dig = "dig",
  antipredator = "vigilance")

#plot
cont1_plot <- fviz_contrib(PCA_results, choice = "var", axes = 1) +
  scale_x_discrete(
    limits = my_order,
    labels = my_labels) +
  labs(title = "Contributors to PC1",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cont2_plot <- fviz_contrib(PCA_results, choice = "var", axes = 2) +
  scale_x_discrete(
    limits = my_order,
    labels = my_labels) +
  labs(title = "Contributors to PC2",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cont3_plot <- fviz_contrib(PCA_results, choice = "var", axes = 3) +
  scale_x_discrete(
    limits = my_order,
    labels = my_labels) +
  labs(title = "Contributors to PC3",
       x = "Variables", y = "Contribution (%)") +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##PC1####

###model PC1####

#create models
model <- lm(PC1 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC1 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC1 ~ cos_cycle*species*treatment , 
                data = scores)

#create reference model
gam_model <- gamm4(
  PC1 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores)

#compare
pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor_matrix <- cbind(pred_gam, pred_m2, pred_m3, pred_m4) #correlation matrix
res <- rcorr(cor_matrix, type = "pearson")

res$r   # correlations
res$P   # p-values
res$n   # n

summary(cos_model)

#remove non-significant interactions
PC1_model <- lm(PC1 ~ cos_cycle+treatment+species, 
                  data = scores)

summary(PC1_model)

#final model
model_PC1.maze <- summary(PC1_model)



###plot PC1####

dat = predict_response(PC1_model, terms=c("cos_cycle", "treatment"))
PC1_plot <- plot(dat, show_data = T, jitter = 0.05) +
  labs(
    title = "Predicted values of PC1",
    x = "Lunar cycle",
    y = "PC1",
    color = "Treatment"
  ) +
  scale_color_manual(
    values = c("#f03b20", "#2c7fb8"),   # Blue and vermillion
    labels = c("Control", "Maze")
  ) +
  scale_fill_manual(
    values = c("#f03b20", "#2c7fb8")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Full", "Waning/Waxing","New"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.spacing.x = unit(3, "lines"))

##PC2####
###model PC2####

#create models
model <- lm(PC2 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC2 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC2 ~ cos_cycle*species*treatment , 
                data = scores)

#create reference model
gam_model <- gamm4(
  PC2 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores
)

#compare
pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor_matrix <- cbind(pred_gam, pred_m2, pred_m3, pred_m4)
res <- rcorr(cor_matrix, 
             type = "pearson")

res$r   # correlations
res$P   # p-values
res$n   # n

summary(cos_model)

#remove non-significant interactions
PC2_model <- lm(PC2 ~ cos_cycle + species + treatment , 
                  data = scores)

#final model
summary(PC2_model)

model_PC2.maze <- summary(PC2_model)


###plot PC2####
dat = predict_response(PC2_model, terms=c("cos_cycle", "treatment"))
PC2_plot <- plot(dat, show_data = T, jitter = 0.05) +
  labs(
    title = "Predicted values of PC2",
    x = "Lunar cycle",
    y = "PC2",
    color = "Treatment"
  ) +
  scale_color_manual(
    values = c("#f03b20", "#2c7fb8"),   # Blue and vermillion
    labels = c("Control", "Maze")
  ) +
  scale_fill_manual(
    values = c("#f03b20", "#2c7fb8")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Full", "Waning/Waxing","New"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.spacing.x = unit(3, "lines"))

##PC3####

###model PC3####

#create models
model <- lm(PC3 ~ sin_cycle*species*treatment +
              cos_cycle*species*treatment , 
            data = scores)
sin_model <- lm(PC3 ~ sin_cycle*species*treatment , 
                data = scores)
cos_model <- lm(PC3 ~ cos_cycle*species*treatment , 
                data = scores)

#create reference model
gam_model <- gamm4(
  PC3 ~
    s(cycle, by = interaction(species,treatment), k = 4) +
    species*treatment,
  data = scores
)

#compare
pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor_matrix <- cbind(pred_gam, pred_m2, pred_m3, pred_m4)
res <- rcorr(cor_matrix, 
             type = "pearson")

res$r   # correlations
res$P   # p-values
res$n   # n

summary(sin_model)

#remove non-significant interactions
PC3_model <- lm(PC3 ~ sin_cycle + species + treatment , 
                  data = scores)

summary(PC3_model)

#final model
model_PC3.maze <- summary(PC3_model)


###plot PC3####

dat = predict_response(PC3_model, terms=c("sin_cycle", "treatment"))
PC3_plot <- plot(dat, show_data = T, jitter = 0.05) +
  labs(
    title = "Predicted values of PC3",
    x = "Lunar cycle",
    y = "PC3",
    color = "Treatment"
  ) +
  scale_color_manual(
    values = c("#f03b20", "#2c7fb8"),   # Blue and vermillion
    labels = c("Control", "Maze")
  ) +
  scale_fill_manual(
    values = c("#f03b20", "#2c7fb8")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Waning", "New/Full", "Waxing"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.spacing.x = unit(3, "lines"))

#time to find the second tray####

##model selection####

#create models
model <- lm(second_tray ~ sin_cycle * species * treatment +
              cos_cycle * species * treatment, 
            data = combined_df)
sin_model <- lm(second_tray ~ sin_cycle * species * treatment, 
                data = combined_df)
cos_model <- lm(second_tray ~ cos_cycle * species * treatment, 
                data = combined_df)

#create reference model
gam_model <- gamm4(
  second_tray ~
    s(cycle, by = interaction(treatment, species), k = 4) +
    treatment * species,
  data = combined_df
)

#compare
pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor_matrix <- cbind(pred_gam, pred_m2, pred_m3, pred_m4)
res <- rcorr(cor_matrix, 
             type = "pearson")

res$r   # correlations
res$P   # p-values
res$n   # n

#log transform
log_model <- lm(log(second_tray) ~ sin_cycle * species * treatment +
                  cos_cycle * species * treatment, 
                data = combined_df)

AIC(model, log_model)

summary(log_model)

#remove non-significant interactions
time_model <- lm(log(second_tray) ~ sin_cycle * treatment +
                   cos_cycle * treatment + species, 
                 data = combined_df)

#final model
summary(time_model)

##get values for maximum and minimum####

# Extract coefficients
b <- coef(time_model)
###control####
a_control <- b["sin_cycle"]
c_control <- b["cos_cycle"]
int_control <- b["(Intercept)"]

# Amplitude
R_control <- sqrt(a_control^2 + c_control^2)

# Angle of maximum (radians)
theta_max_control <- atan2(a_control, c_control)
if(theta_max_control < 0) theta_max_control <- theta_max_control + 2*pi

# Angle of minimum
theta_min_control <- theta_max_control + pi
if(theta_min_control > 2*pi) theta_min_control <- theta_min_control - 2*pi

# Convert to lunar days
cycle_length <- 29.53

day_max_control <- theta_max_control/(2*pi) * cycle_length
day_min_control <- theta_min_control/(2*pi) * cycle_length

# Predicted times (minutes)
max_control <- exp(int_control + R_control)
min_control <- exp(int_control - R_control)



###maze####
a_maze <- b["sin_cycle"] + b["sin_cycle:treatmentmaze"]
c_maze <- b["cos_cycle"] + b["treatmentmaze:cos_cycle"]
int_maze <- b["(Intercept)"] + b["treatmentmaze"]

R_maze <- sqrt(a_maze^2 + c_maze^2)

theta_max_maze <- atan2(a_maze, c_maze)
if(theta_max_maze < 0) theta_max_maze <- theta_max_maze + 2*pi

theta_min_maze <- theta_max_maze + pi
if(theta_min_maze > 2*pi) theta_min_maze <- theta_min_maze - 2*pi

day_max_maze <- theta_max_maze/(2*pi) * cycle_length
day_min_maze <- theta_min_maze/(2*pi) * cycle_length

max_maze <- exp(int_maze + R_maze)
min_maze <- exp(int_maze - R_maze)
###summary####
results <- data.frame(
  treatment = c("Control", "Maze"),
  min_minutes = c(min_control, min_maze),
  min_cycle_day = c(day_min_control, day_min_maze),
  max_minutes = c(max_control, max_maze),
  max_cycle_day = c(day_max_control, day_max_maze)
)

print(results, digits = 3)
##plot####

dat = predict_response(time_model, terms=c("cos_cycle","treatment"))

time_plot.cos <- plot(dat, show_data = TRUE, jitter = 0.05) +
  labs(
    title = "Predicted time to second tray (min)",
    x = "Lunar cycle",
    y = "Time to second tray (min)",
    color = "Treatment"
  ) +
  scale_color_manual(
    values = c("#f03b20", "#2c7fb8"),   # Blue and vermillion
    labels = c("Control", "Maze")
  ) +
  scale_fill_manual(
    values = c("#f03b20", "#2c7fb8")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Full", "Waning/Waxing","New"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.spacing.x = unit(3, "lines"))


dat = predict_response(time_model, terms=c("sin_cycle","treatment"))
time_plot.sin <- plot(dat, show_data = TRUE, jitter = 0.05) +
  labs(
    title = "Predicted time to second tray (min)",
    x = "Lunar cycle",
    y = "Time to second tray (min)",
    color = "Treatment"
  ) +
  scale_color_manual(
    values = c("#f03b20", "#2c7fb8"),   # Blue and vermillion
    labels = c("Control", "Maze")
  ) +
  scale_fill_manual(
    values = c("#f03b20", "#2c7fb8")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Waning", "New/Full", "Waxing"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(3, "lines"))

#seed intake####

##create data set####
gud_df <- pivot_longer(combined_df,
                       cols = c(seeds_maze,seeds_control), 
                       names_to = "tray",          
                       values_to = "gram")
gud_df <- gud_df[which(gud_df$trial_length == 60),]
gud_df$treatment <- as.factor(gud_df$treatment)
gud_df$tray <- as.factor(gud_df$tray)
gud_df$species <- as.factor(gud_df$species)

##model selection####

#create models
model <- lm(gram ~ sin_cycle * species * treatment * tray +
              cos_cycle * species * treatment * tray, 
            data = gud_df)
sin_model <- lm(gram ~ sin_cycle * species * treatment * tray, 
                data = gud_df)
cos_model <- lm(gram ~ cos_cycle * species * treatment * tray, 
                data = gud_df)

#create reference model
gam_model <- gamm4(
  gram ~
    s(cycle, by = interaction(treatment, species, tray), k = 4) +
    treatment * species * tray,
  data = gud_df
)

#compare
pred_gam  <- predict(gam_model$gam, type = "response")
pred_m2   <- predict(model, type = "response")
pred_m3   <- predict(sin_model, type = "response")
pred_m4   <- predict(cos_model, type = "response")

cor_matrix <- cbind(pred_gam, pred_m2, pred_m3, pred_m4)
res <- rcorr(cor_matrix, 
             type = "pearson")

res$r   # correlations
res$P   # p-values
res$n   # n

AIC(sin_model,cos_model,model)

summary(sin_model)

#remove non-significant interactions
seed_model <- lm(gram ~ sin_cycle*tray + tray*treatment + species, 
            data = gud_df)

#final model
summary(seed_model)

#effects of tray while considering treatment
emmeans(seed_model, ~ tray*treatment)

#effects of lunar cycle on the two trays
emtrends(seed_model, ~tray, var = "sin_cycle") 

##plot####
dat = predict_response(seed_model, terms=c("sin_cycle", "tray", "treatment"))
dat$group <- fct_relevel(dat$group, "seeds_maze", "seeds_control")

seed_plot <- plot(dat, show_data = TRUE, jitter = 0.05) +
  labs(
    title = "Predicted seed intake (gram)",
    x = "Lunar cycle",
    y = "Seeds foraged (g)",
    color = "Food tray"
  ) +
  scale_color_manual(
    values = c("#D55E00", "#0072B2"),   # Blue and vermillion
    labels = c("Second tray", "First tray")
  ) +
  scale_fill_manual(
    values = c("#2c7fb8", "#D55E00")) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Waning", "New/full", "Waxing"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(3, "lines"))

#save models & plots####

write.csv(PCA_loadings,'PCA_loadings.csv')
write.csv(PCA_eigenvalues,'PCA_eigenvalues.csv')
write.csv(model_PC1.maze$coefficients, 'model_PC1.csv')
write.csv(model_PC2.maze$coefficients, 'model_PC2.csv')
write.csv(model_PC3.maze$coefficients, 'model_PC3.csv')
write.csv(sum.time_model$coefficients, 'time_model.csv')
write.csv(seed_model$coefficients, 'seed.intake_model.csv')

ggsave(filename = "cont1_plot.png",
  plot = cont1_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "cont2_plot.png",
  plot = cont2_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "cont3_plot.png",
  plot = cont3_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "PC1_plot.png",
  plot = PC1_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "PC2_plot.png",
  plot = PC2_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "PC3_plot.png",
  plot = PC3_plot,
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "seed_plot.png",
  plot = seed_plot,
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "time_plot.cos.png",
  plot = time_plot.cos,
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600)

ggsave(filename = "time_plot.sin.png",
  plot = time_plot.sin,
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600)