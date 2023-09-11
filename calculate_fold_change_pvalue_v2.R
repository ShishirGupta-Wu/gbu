

setwd("C:/Users/Shishir/Desktop/Shared@H/HHU_PROJECTS/PR2-Iron_LEAFAREA/4-Results/stats_revised/Anova_test/t-test-error")

library(dplyr)
library(readr)
library(janitor)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(purrr)

leaf_expression <- read_csv("Leaf_expression.csv", col_names = TRUE) 
#View(Leaf_expression)

# Remove rows where Day is equal to "Day 0"
leaf_expression <- leaf_expression[leaf_expression$Day != "Day 0", ]
colnames(leaf_expression)

# Convert variables to factors
leaf_expression$Treatment <- as.factor(leaf_expression$Treatment)
leaf_expression$Line <- as.factor(leaf_expression$Line)
leaf_expression$Day <- as.factor(leaf_expression$Day)

glimpse(leaf_expression)

#log2 value +1 transformation
#leaf_expression <- leaf_expression %>%   mutate_at(vars(-Datapoints, -Line, -Day, -Treatment), ~ log2(. + 1))

####################
#log2 transformation
#####################

leaf_expression <- leaf_expression %>%   mutate_at(vars(-Datapoints, -Line, -Day, -Treatment), ~ log2(.))



#######################
# Calculate fold change
########################

calculate_fold_change <- function(data, variable, reference_treatment, control) {
  data %>%
    group_by(Line, Day) %>%
    summarize(
      fold_change = mean({{ variable }}[Treatment == reference_treatment]) / mean({{ variable }}[Treatment == control], na.rm = TRUE)
    )
}


gene_names <- leaf_expression %>%  select(-Line, -Datapoints, -Treatment, -Day) %>%  colnames()
gene_names

gene_names_array <- as.array(gene_names)
gene_names_array

fold_changes <- bind_rows(lapply(gene_names_array, function(variable) {
  calculate_fold_change(data = leaf_expression, variable = !!sym(variable), reference_treatment = "Fe0", control = "Fe")
}), .id = "variable") %>%
  pivot_wider(names_from = "variable", values_from = "fold_change")

colnames(fold_changes)[-c(1:2)] <- gene_names_array

fold_changes

print(fold_changes, n = 20)




#######################
# Adding p-values
########################



gene_names <- leaf_expression %>%  select(-Line, -Datapoints, -Treatment, -Day) %>%  colnames()
gene_names

gene_names_array <- as.array(gene_names)
gene_names_array


calculate_fold_change <- function(data, variable, reference_treatment, control) {
  data %>%
    group_by(Line, Day) %>%
    summarize(
      fold_change = mean(.data[[variable]][Treatment == reference_treatment]) / mean(.data[[variable]][Treatment == control], na.rm = TRUE)
    )
}




fold_changes <- bind_rows(lapply(gene_names, function(variable) {
  result <- calculate_fold_change(data = leaf_expression, variable = variable, reference_treatment = "Fe0", control = "Fe")
  if (length(unique(na.omit(result$fold_change))) <= 1) {
    result$p_value <- NA
  } else {
    tryCatch({
      t_test_result <- t.test(result$fold_change, na.rm = TRUE)
      result$p_value <- t_test_result$p.value
    }, error = function(e) {
      result$p_value <- NA
    })
  }
  result
}), .id = "variable") %>%
  pivot_wider(names_from = "variable", values_from = c("fold_change", "p_value"))


colnames(fold_changes)

colnames(fold_changes)[-c(1:2)] <- c(gene_names, paste0(gene_names, "_pval"))
colnames(fold_changes)

print(fold_changes, n = 20)


###############
#Adding p-value stars
###############





# Define the p-value thresholds and corresponding stars
p_value_thresholds <- c(0.001, 0.01, 0.05)
p_value_stars <- c("***", "**", "*")

# Add p-value stars column
fold_changes <- fold_changes %>%
  mutate(across(ends_with("_pval"), ~ case_when(
    . <= p_value_thresholds[1] ~ p_value_stars[1],
    . <= p_value_thresholds[2] ~ p_value_stars[2],
    . <= p_value_thresholds[3] ~ p_value_stars[3],
    TRUE ~ ""
  ), .names = "{.col}_stars"))

# Print the modified fold_changes
print(fold_changes, n = Inf)

# Open the modified fold_changes in a new viewer window
View(fold_changes)














