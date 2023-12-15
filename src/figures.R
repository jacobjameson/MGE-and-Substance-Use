#-------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate figures
#-------------------------------------------------------------------------

# Create dataset ----------------------------------------------------------

source('src/Construct Analytical Dataset.R')

# load packages ----------------------------------------------------------
library(marginaleffects)
library(ggeffects)
library(ggsci)
library(ggnetwork)
library(intergraph)
library(igraph)
library(network)
library(GGally)


##########################################################################
# FIGURE 1
#
# Figure 1: Social network map of friendship nominations for 
# 1 of 142 schools included in analysis
##########################################################################

data_path <- '~/Add Health Data'
inschool_path <-  paste0(data_path, 
                         '/Wave I In-School Questionnaire Data')

inschool <- read_xpt(paste0(inschool_path, '/Inschool.xpt'))
inschool <- inschool[,c('AID', 'SQID', 'SSCHLCDE')]

inschool <- inschool %>%
  filter(SQID != '', AID != '')

friend.df <- read_xpt(paste0(data_path,
                             '/Wave I In-School Friendship Nominations/sfriend.xpt'))

friend.df <- friend.df %>%
  filter(SQID != '999999')

friend.df <- merge(friend.df, inschool, by='SQID')
friend.df[] <- lapply(friend.df, as.character)
names(friend.df) <- tolower(names(friend.df))

friend.vars <- c('mf1aid', 'mf2aid', 'mf3aid', 'mf4aid', 'mf5aid', 
                 'ff1aid', 'ff2aid', 'ff3aid', 'ff4aid', 'ff5aid')

id.replace <- c('77777777', '99999999', '88888888', '99959995')


friend.df <- friend.df %>%
  mutate(mf1aid = ifelse(mf1aid %in% id.replace, NA, mf1aid),
         mf2aid = ifelse(mf2aid %in% id.replace, NA, mf2aid),
         mf3aid = ifelse(mf3aid %in% id.replace, NA, mf3aid),
         mf4aid = ifelse(mf4aid %in% id.replace, NA, mf4aid),
         mf5aid = ifelse(mf5aid %in% id.replace, NA, mf5aid),
         ff1aid = ifelse(ff1aid %in% id.replace, NA, ff1aid),
         ff2aid = ifelse(ff2aid %in% id.replace, NA, ff2aid),
         ff3aid = ifelse(ff3aid %in% id.replace, NA, ff3aid),
         ff4aid = ifelse(ff4aid %in% id.replace, NA, ff4aid),
         ff5aid = ifelse(ff5aid %in% id.replace, NA, ff5aid))


friend.long <- friend.df[,c('aid', 'sschlcde', 'mf1aid', 'mf2aid', 'mf3aid', 'mf4aid',
                              'mf5aid', 'ff1aid','ff2aid', 'ff3aid', 'ff4aid',
                              'ff5aid')]  %>% 
  pivot_longer(cols = c(mf1aid, mf2aid, mf3aid, mf4aid, mf5aid, 
                        ff1aid, ff2aid, ff3aid, ff4aid, ff5aid), 
               names_to = "friend_type", values_to = "friend_aid") %>% 
  filter(!is.na(friend_aid))


network <- merge(
  friend.long, 
  analytical_dataset %>% 
    mutate(high.ge = ifelse(w1.GE_male_std >= 1, 'Males with High Adolescent GE', 'Other Students')) %>%
  filter(sschlcde == '077', is.na(w1.GE_male_std)==F) %>%
    dplyr::select(aid, in_sample, w1.GE_male_std, w4.GE_male_std, delta_w1_w4_GE, high.ge, sschlcde),
  by=c('aid'), all.x = T)

network <- network[network$aid != network$friend_aid, ]
network <- network[network$friend_aid %in% network$aid, ]

network <- filter(network, is.na(w1.GE_male_std) == F)

net <- list(nodes=network[c('aid', 'friend_aid', 'high.ge')], 
            edges=network[c('high.ge', 'aid', 'friend_aid')])

# create node attribute data
net.cet <- as.character(net$nodes$high.ge)
names(net.cet) = net$nodes$aid
edges <- net$edges

# create network
net.net <- edges[, c("aid", "friend_aid") ]
net.net <- network(net.net, directed = TRUE)

# create GE type node attribute
net.net %v% "high.ge" <- net.cet[network.vertex.names(net.net) ]

set.seed(10312016)
ggnet2(net.net, color = "high.ge", 
       palette = c('Males with High Adolescent GE' = "red", 'Other Students' = 'grey50'), size = 'indegree',
       arrow.size = 3, arrow.gap = 0.02, alpha = .9, edge.size = 0.5,
       edge.alpha = 0.9,  mode = "fruchtermanreingold",  
       color.legend = "") +
  theme_bw() + theme_blank()  + 
  theme(legend.position = "bottom", text = element_text(size = 20),
        plot.caption = element_text(hjust = 0)) + guides(size=F) + 
  labs(title='Social Network Mapping of Friendship Nominations for 1 of 142 Schools Included in the Analysis\n', 
       caption = str_wrap("\nThe size of the dots represent the number of friendship nominations by others. 
                          The red dots represent males with gender expression scores that are more than one 
                          standard deviation greater than that of the male average. Directed edges are used 
                          to indicate the direction of friendship nominations. Edges that are bi-directional 
                          indicated friendship reciprocity.\n\n This particular 
                          network represents the frienship nominations of 624 students from the Wave I In-School 
                          Friendship Nominations data.", 150)) +
  guides(color=guide_legend(keyheight=0.5,default.unit="inch",override.aes = list(size=10)))

ggsave("outputs/figures/Figure 1.png", width = 15, height = 10, bg = 'white')


##########################################################################
# FIGURE 2
#
# Figure 2: Marginal effects of gender expression (GE) in adolescence 
# and adulthood on predicted young adult substance use
##########################################################################

# Define a function to fit the model
fit_model <- function(outcome, predictor, data) {
  formula <- as.formula(paste(outcome, "~", predictor, 
                              "+ race + pseudo.gpa + sespc_al + nhood1_d"))
  glm(formula, data = subset(data, in_sample == 1), 
      weights = weights, family = 'quasibinomial')
}

# Define a function to get predictions
get_predictions <- function(model, predictor) {
  ggpredict(model, terms = c(paste(predictor, "[-2:2]")), vcov.type = "HC0", 
            vcov.args = list(cluster = subset(analytical_dataset, 
                                              in_sample == 1)$cluster)) %>%
    as.data.frame() %>%
    dplyr::select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high) %>%
    filter(!is.na(xvals))
}

# Define outcomes and their corresponding predictors
outcomes <- c("w4.cigarettes.bin.30", "w4.marijuana.bin.30", 
              "w4.drunk.bin.30", "w4.fav.bin.30", "w4.prescription")

predictors <- c("w4.GE_male_std", "w4.GE_male_std", "w4.GE_male_std", 
                "w4.GE_male_std", "w4.GE_male_std")

labels <- c("Cigarette Smoking", "Marijuana Use", "Excessive Alcohol Use",
            "Recreational Drug Use", "Prescription Drug Misuse")

data_list <- list()

# Loop through outcomes and predictors to fit models and get predictions
for (i in 1:length(outcomes)) {
  model <- fit_model(outcomes[i], predictors[i], analytical_dataset)
  df <- get_predictions(model, predictors[i])
  df$Variables <- labels[i]
  data_list[[i]] <- df
}

# Combine data frames
data <- do.call(rbind, data_list)
data$group <- 'Model 2: Young Adult (Wave IV) MGE'

# Repeat for the second set of predictors
predictors <- c("w1.GE_male_std", "w1.GE_male_std", "w1.GE_male_std", 
                "w1.GE_male_std", "w1.GE_male_std")

data_list_2 <- list()

for (i in 1:length(outcomes)) {
  model <- fit_model(outcomes[i], predictors[i], analytical_dataset)
  df <- get_predictions(model, predictors[i])
  df$Variables <- labels[i]
  data_list_2[[i]] <- df
}

data_2 <- do.call(rbind, data_list_2)
data_2$group <- 'Model I: Adolsecent (Wave I) MGE'

# Combine both datasets
data <- rbind(data_2, data)

# Factorize variables for ordering
data$Variables <- factor(data$Variables, levels = c("Excessive Alcohol Use", "Cigarette Smoking", 
                                                    'Prescription Drug Misuse', 'Marijuana Use', 
                                                    'Recreational Drug Use'))

data$group <- factor(data$group, levels = c('Model I: Adolsecent (Wave I) MGE', 
                                            'Model 2: Young Adult (Wave IV) MGE'))


ggplot(data, aes(x=xvals,fill=Variables)) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5) +
  geom_line(aes(y=coef)) + 
  scale_fill_nejm() +
  facet_wrap(~group, ncol = 2) +
  theme_bw() +
  geom_vline(xintercept =0, color='darkred') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits =c(0, 0.9)) +
  labs(fill = "Wave IV Behavior",
       y = "Predicted Probability\n", 
       x='\nStandardized Gender Expression\n',
       title = 'Marginal effects of male gender expression (MGE) in adolescence and young adulthood on predicted adult substance use\n',
       caption = str_wrap('\n\nIncreases in adolescent  (Wave I) MGE and young adult (Wave IV) MGE 
                          are associated with higher predicted probabilities of young adult 
                          substance use. Note: Adolescence, Participants aged 12-18; Young 
                          adulthood, Participants aged 24-32', 140)) +
  theme(axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.x = element_text(color = 'black',size = 18),
        axis.title.y = element_text(color = 'black', size = 18),
        plot.title = element_text(color = "black", size = 18, hjust = 0),
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        legend.position = 'right',
        legend.title = element_text(size=18),
        strip.text.x = element_text(color = 'black', size = 18, face = "bold"),
        legend.text = element_text(color = "black", size = 18))


ggsave("outputs/figures/Figure 2.png", width = 15, height = 9, bg = 'white')
ggsave("outputs/figures/Figure 2.pdf", width = 15, height = 9)

##########################################################################

##########################################################################
# FIGURE 3
#
# Figure 3: Marginal effects of adolescent-to-adult change in
# gender expression (GE) on predicted adult substance use
##########################################################################

# Define a function to fit the model
fit_model <- function(outcome, predictor, data) {
  formula <- as.formula(paste(outcome, "~ delta_w1_w4_GE + race + pseudo.gpa + sespc_al + nhood1_d", 
                              ifelse(outcome != "w4.prescription", paste("+", predictor), "")))
  glm(formula, data = subset(data, in_sample == 1), weights = weights, family = 'quasibinomial')
}

# Define a function to get predictions
get_predictions <- function(model) {
  ggpredict(model, terms = c('delta_w1_w4_GE [-2:2]'), vcov.type = "HC0", 
            vcov.args = list(cluster = subset(final.df, in_sample == 1)$cluster)) %>%
    as.data.frame() %>%
    dplyr::select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high) %>%
    filter(!is.na(xvals))
}

# Define outcomes and their corresponding predictors
outcomes <- c("w4.cigarettes.bin.30", "w4.marijuana.bin.30", "w4.drunk.bin.30", "w4.fav.bin.30", "w4.prescription")
predictors <- c("w1.cigarettes", "w1.marijuana", "w1.drunk", "w1.recreational", NULL)
labels <- c("Cigarette Smoking", "Marijuana Use", "Excessive Alcohol Use", "Recreational Drug Use", "Prescription Drug Misuse")

data_list <- list()

# Loop through outcomes and predictors to fit models and get predictions
for (i in 1:length(outcomes)) {
  model <- fit_model(outcomes[i], predictors[i], analytical_dataset)
  df <- get_predictions(model)
  df$Variables <- labels[i]
  data_list[[i]] <- df
}

# Combine data frames
data <- do.call(rbind, data_list)
data$group <- 'Model 3: Change in Adolescent to Adult GE'


# Factorize variables for ordering
data$Variables <- factor(data$Variables, levels = c("Excessive Alcohol Use", "Cigarette Smoking", 
                                                    'Prescription Drug Misuse', 'Marijuana Use', 
                                                    'Recreational Drug Use'))

ggplot(data, aes(x=xvals,fill=Variables)) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5) +
  geom_line(aes(y=coef)) + scale_fill_nejm() +
  geom_vline(xintercept =0, color='darkred') +
  annotate("text", x = -1, y = 0.82, color ='black', fontface =2,
           label = "∆ MGE is Decreasing\n←", size = 6) +
  annotate("text", x = 1, y = 0.82, color ='black', fontface =2,
           label = "∆ MGE is Increasing\n→", size = 6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits =c(0.02, 0.85)) +
  theme_minimal() +
  labs(fill = "Wave IV Behavior",
       y = "Predicted Probability\n", 
       x='\nChange in Standardized Male Gender Expression \nfrom Adolescence to Young Adulthood\n',
       title = str_wrap('Marginal effects of adolescent-to-young-adult change in male gender expression (MGE) on predicted young adult substance use\n\n', 100),
       caption = str_wrap('\n\nIncreases in relative MGE are associated with higher predicted 
                          probabilities of young adult substance use behavior. Note: 
                          Adolescence, Participants aged 12-18; Young adulthood, Participants aged 24-32', 120)) +
  theme(axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.x = element_text(color = 'black',size = 18),
        axis.title.y = element_text(color = 'black', size = 18),
        plot.title = element_text(color = "black", size = 18, hjust = 0),
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        legend.position = 'right',
        legend.title = element_text(size=18),
        strip.text.x = element_text(color = 'black', size = 18, face = "bold"),
        legend.text = element_text(color = "black", size = 18))


ggsave("outputs/figures/Figure 3.png", width = 14, height = 10, bg = 'white')

##########################################################################

##########################################################################
# Appendix 9
#
# Appendix 9. Incidence rate-ratios from negative biomaial regression 
# models assessing impact of gender expression over time on adult
# substance use over prior 30-days
##########################################################################

# Define the models
models <- list(
  w4.cigarettes.cont.30 = "w1.cigarettes",
  w4.marijuana.cont.30 = "w1.marijuana",
  w4.drunk.cont.30 = "w1.drunk",
  w4.fav.cont.30 = "w1.recreational")

# Define the list of predictors
predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
base_predictors <- "race + pseudo.gpa + sespc_al + nhood1_d"

nbregs <- data.frame()

# Loop through each predictor
for (predictor in predictors) {
  
  # Loop through each model and run regression for current predictor
  for (var in names(models)) {
    
    # Construct the formula
    unique_predictor <- models[[var]]

    formula_str <- paste(var, "~", base_predictors, "+", predictor, "+", unique_predictor)
    formula_obj <- as.formula(formula_str)
    
    # Run regression
    mod <- glm.nb(formula = formula_obj, data = subset(analytical_dataset, in_sample == 1),
                  weights = weights)

    robust <- coeftest(x = mod, vcov = vcovCL(mod, type = "HC0", cluster = ~ cluster))
    ci  <- confint(robust)
    
    robust <- robust %>% 
      tidy() %>% 
      as.data.frame() %>%
      mutate(lwr = as.data.frame(ci)[,1],
             upr = as.data.frame(ci)[,2],
             outcome = var) %>%
      filter(term %in% predictors)
    
    # Get IRR
    robust$estimate <- exp(robust$estimate)
    robust$std.error <- exp(robust$std.error)
    robust$lwr <- exp(robust$lwr)
    robust$upr <- exp(robust$upr)
    
    nbregs <- bind_rows(nbregs, robust)
  }
}


nbregs %>%
  mutate(vars = factor(case_when(
    term == 'w1.GE_male_std' ~ 'Wave I GE Measure',
    term == 'w4.GE_male_std' ~ 'Wave IV GE Measure',
    term == 'delta_w1_w4_GE' ~ 'Change in Wave I to Wave IV GE'),
    levels = c('Change in Wave I to Wave IV GE', 'Wave IV GE Measure', 'Wave I GE Measure')),
         outcomes = factor(case_when(
           outcome == 'w4.cigarettes.cont.30' ~ 'Wave IV 30 Day Cigarette Use',
           outcome == 'w4.marijuana.cont.30' ~ 'Wave IV 30 Day Marijuana Use',
           outcome == 'w4.drunk.cont.30' ~ 'Wave IV 30 Day Alcohol Use',
           outcome == 'w4.fav.cont.30' ~ 'Wave IV 30 Day Recreational Drug Use'),
    levels = c('Wave IV 30 Day Cigarette Use', 'Wave IV 30 Day Marijuana Use',
               'Wave IV 30 Day Alcohol Use', 'Wave IV 30 Day Recreational Drug Use'))) %>%
  mutate(zero = ifelse(lwr<1 & upr>1, "p-value > 0.05", "p-value < 0.05" )) %>% 
  ggplot() + 
  geom_linerange(mapping = aes(x = vars, ymin = lwr, ymax = upr), size = 0.75) +
  geom_point(mapping = aes(x = vars, y = estimate, color = zero), size = 3.5) +
  geom_hline(mapping = aes(yintercept = 1), linetype = "dashed", color = 'darkred') + 
  facet_wrap(~outcomes, nrow=4) + coord_flip() + 
  labs(x = "", 
       y = "\nIncidence Rate Ratio from Negative Binomial Model", 
       title = 'Incidence rate-ratios from negative binomial regression models assessing impact of male gender \nexpression (MGE) over time on young adult substance use over prior 30-days\n',
       color = "Significance") + 
  scale_color_d3() +  theme_bw() +
  theme(strip.text.x = element_text(face ="bold", size =18),
        axis.text = element_text(colour = "black", size=18),
        legend.position = 'bottom',
        axis.title = element_text(size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_text (colour ="black", size =18 ,
                                     face ="bold"),
        legend.background = element_rect (size =0.5 , linetype ="solid",
                                          colour ="black"),
        plot.title = element_text(size =18))

ggsave("outputs/figures/Appendix 9.png", width = 16, height = 10, bg = 'white')

