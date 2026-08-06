## Create contact matrix using EpiStorm-Mix
## Code modified from: https://github.com/epistorm/Epistorm-Mix/
## Accessed 23 Dec 2025

## Load & merge data -------------

#respondent-level data
resp.url <- "https://raw.githubusercontent.com/epistorm/Epistorm-Mix/refs/heads/main/tutorials/data-analysis/data/imputed_resp.csv" 
resp_imp <- read.csv(resp.url)

# contact-level data
cont.url <- "https://raw.githubusercontent.com/epistorm/Epistorm-Mix/refs/heads/main/tutorials/data-analysis/data/imputed_cont.csv" 
cont_imp <- read.csv(cont.url)

# merge
merged.df <- left_join(resp_imp, cont_imp, by = "resp_id")


## Contact matrix with 5 age groups -----------------

# Set age brackets and labels
age.brks <- c(0,5, 18, 50, 65, Inf)
grp.lbls <- c("0-4", "5-17", "18-49", "50-64","65+")

merged.df <- merged.df %>% 
              mutate(resp_agegrp = cut(resp_age, breaks = age.brks, labels = grp.lbls, right = FALSE),
                     cont_agegrp = cut(cont_age, breaks = age.brks, labels = grp.lbls, right = FALSE))

resp_imp <- resp_imp %>% mutate(resp_agegrp = cut(resp_age, breaks = age.brks, labels = grp.lbls, right = FALSE))

# weighted N of resp
age_matrix_n <- resp_imp %>% group_by(resp_agegrp) %>%
                  summarise(n = sum(weight, na.rm = TRUE))

# all possible age group combinations
nd <- expand.grid(resp_agegrp = levels(merged.df$resp_agegrp),
                  cont_agegrp = levels(merged.df$cont_agegrp))
predict.m <- as.data.frame(nd)

# obtain the weighted sum of contacts per respondent-contact group
age_matrix <- merged.df %>% 
                filter(!is.na(cont_agegrp)) %>%  # exclude NAs (contacts with missing ages)
                arrange(cont_agegrp) %>%
                group_by(cont_agegrp, resp_agegrp) %>%
                # each row represents a contact, so just sum up weight (x 1) for each row
                summarise(weighted_sum = sum(weight, na.rm = T)) 

# combine weighted sum with the data for predictions 
# (would have NAs if there were any missing combinations)
age_matrix <- left_join(predict.m, age_matrix)

# order (respondent first)
age_matrix <- age_matrix[order(age_matrix[,"resp_agegrp"], age_matrix[,"cont_agegrp"]),]

# combine with total (weighted) # respondents
age_matrix <- left_join(age_matrix, age_matrix_n)

# average number of contacts for the matrix: divide respondent-contact pair by total # respondents
age_matrix$ncont <- age_matrix$weighted_sum/age_matrix$n   

# convert to matrix format
## transform dataframe into matrix form
contact_mat <- age_matrix %>% 
                dplyr::select("resp_agegrp", "cont_agegrp", "ncont") %>% 
                pivot_wider(names_from = "cont_agegrp", values_from = "ncont") %>%
                select(-resp_agegrp)

                     