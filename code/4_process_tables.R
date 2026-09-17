xfun::pkg_attach(c("tidyverse", "gt", "here", "knitr", "kableExtra"))

tabpath <- "tables"
filepath <- ""

source(here(filepath, 'code/0_functions.R'))

# scenario information
scenario_tab <- read.csv(here(filepath,'inputs/scenarios.csv'), row.names = 1)
epi_scenarios <- names(scenario_tab)

scenario_short <- unique(str_split_i(epi_scenarios, pattern = "_", i = 3))
scenario_labels <- c("Baseline", "No transmission impact", "Small transmission impact", "Large transmission impact")

age_labels <- c("0-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs")


# Output cleaning ------

# load results from main scenario analysis (1_run_model_original.R)
load(here(filepath, paste0("out/original_100.RData")) )

tmp_ar <- out %>% 
  filter(time == max(time)) %>%
  gather(var, val, contains("H")) %>%
  mutate(agegrp = substr(var, nchar(var), nchar(var))) %>%
  select(sim, scenario, var, agegrp, val) %>% 
  group_by(sim, scenario, agegrp) %>% 
  summarize(val = sum(val)) %>% ungroup() %>%
  separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
  select(-X)

df0 <- tmp_ar %>% spread(scenario, val) %>%
  # calculate absolute difference and percentage difference  
  mutate(across(.cols = scenario_short, 
                .fns = list(total   = function(x) baseline - x, 
                            indirect = function(x) notrans - x),
                .names = "{.col}_{.fn}")) 

df <- df0 %>% select(sim, severity, agegrp, baseline, contains("total"), -baseline_total)


## Total effects by age group ----------

# Compute percentages within each simulation, then average ---
df_pct <- df %>%
          group_by(sim, severity) %>%
          mutate(across(.cols = contains("_total"), 
                        .fns = list(pct_all   = function(x)  x/sum(x) *100,
                                    pct_within = function(x) x / baseline * 100),
                        .names = "{.col}_{.fn}")) %>%
          ungroup() %>%
          group_by(severity, agegrp) %>%
          summarise(across(ends_with("_pct_all") | ends_with("_pct_within"), mean),
                    .groups = "drop")  %>%
          arrange(severity, agegrp) 

# Compute totals column for % prevented within group ---
# = mean across sims of (sum of prevented / sum of baseline)
df_totals <- df %>%
  group_by(sim, severity) %>%
  summarise(across(contains("_total"),
                   .fns = list(pct_within = function(x) sum(x) / sum(baseline) * 100),
                   .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  group_by(severity) %>%
  summarise(across(ends_with("_pct_within"), mean), .groups = "drop")

# Append total rows ---
total_rows <- df_totals %>%
    mutate(
      agegrp                   = 6,
      notrans_total_pct_all    = 100,
      smalltrans_total_pct_all = 100,
      bigtrans_total_pct_all   = 100
    ) %>%
    select(severity, agegrp,
           notrans_total_pct_all, notrans_total_pct_within,
           smalltrans_total_pct_all, smalltrans_total_pct_within,
           bigtrans_total_pct_all, bigtrans_total_pct_within)
  
tbl <- bind_rows(df_pct, total_rows %>% mutate(agegrp = as.character(agegrp))) %>%
  arrange(severity, agegrp) %>%
  mutate(agegrp_txt = factor(c(age_labels, "Total")[as.numeric(agegrp)], levels = c(age_labels, "Total")),
         severity = paste(str_to_title(severity), "severity season"))

  
#  Build initial table ---

tab_age <- 
  tbl %>%
  select(severity, agegrp_txt,
         notrans_total_pct_within, notrans_total_pct_all,
         smalltrans_total_pct_within, smalltrans_total_pct_all,
         bigtrans_total_pct_within, bigtrans_total_pct_all) %>%
  gt(groupname_col = "severity", rowname_col = "agegrp_txt") %>%
  tab_spanner(label = "No transmission impact",    columns = c(notrans_total_pct_within,    notrans_total_pct_all))    %>%
  tab_spanner(label = "Small transmission impact", columns = c(smalltrans_total_pct_within, smalltrans_total_pct_all)) %>%
  tab_spanner(label = "Large transmission impact",   columns = c(bigtrans_total_pct_within,   bigtrans_total_pct_all))   %>%
  cols_label(
    notrans_total_pct_within      = "% prevented within group",
    notrans_total_pct_all         = "% of all prevented",
    smalltrans_total_pct_within   = "% prevented within group",
    smalltrans_total_pct_all      = "% of all prevented",
    bigtrans_total_pct_within     = "% prevented within group",
    bigtrans_total_pct_all        = "% of all prevented"
  ) %>%
    fmt_number(columns = everything(), decimals = 1, suffix = "%") %>%
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_stub(rows = agegrp_txt == "Total")
    ) %>%
    tab_footnote("% of all prevented: number prevented in age group / total number prevented.") %>%
    tab_footnote("% prevented within group: number prevented in age group / baseline number in age group.") %>%
    opt_row_striping() 


# Re-format ---

tab_age %>% 
  as_latex() %>%
  as.character() %>%
  
  # column line formatting
  gsub("^ & \\\\multicolumn", "\\\\multicolumn", ., perl = TRUE) %>%

  # reformat column headers and text size
  gsub("\\\\% prevented within group", "\\\\makecell{\\\\% prevented \\\\\\\\ within group}", .) %>%
  gsub("\\\\% of all prevented", "\\\\makecell{\\\\% of all \\\\\\\\ prevented}", .) %>%
  gsub("\\\\fontsize\\{.*?\\}\\{.*?\\}\\\\selectfont", "", .) %>%
  
  # reformat footnote
  gsub("\\\\begin\\{minipage\\}\\{\\\\linewidth\\}", 
       "\\\\begin{minipage}{\\\\linewidth}\n\\\\footnotesize{", .) %>%
  gsub("\\\\end\\{minipage\\}", "}\n\\\\end{minipage}", .) %>%
  
  # other text formatting
  gsub("65\\+", "$\\\\geq$65", .) %>%
  
  # add row spacing
  gsub("(\\\\begin\\{tabular\\*\\})", "\\\\renewcommand{\\\\arraystretch}{1.4}\n\\1", .) %>%
  
  # add line above Total rows
  gsub("(\\{\\\\bfseries Total\\})", "\\\\midrule\n\\1", .) %>%
  
  writeLines(here(tabpath, "prevented_hospitalizations.tex"))


## Total undergoing prompt treatment ---------

out_all <- combine_out(out) %>% 
  select(time, scenario, sim, contains("all")) %>% 
  separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels))


tx_tab <- 
  out_all %>%
  filter(time == max(time), scenario != "baseline") %>% 
  select(time, severity, scenario, sim, count_Ia1_all, count_Ia2_all) %>%
  mutate(tx = count_Ia1_all + count_Ia2_all) %>%
  group_by(severity, scenario) %>%
  summarize(mean = mean(tx), 
            lcl = quantile(tx, 0.025),
            ucl = quantile(tx, 0.975)) %>%
  ungroup() %>%
  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels),
         severity = paste(str_to_title(severity), "severity season")) %>%
  mutate(across(c(mean, lcl, ucl), function(x) prettyNum(signif(x, 3), big.mark = ","))) %>% 
  mutate(estimate = paste0(mean, " (", lcl, "--", ucl, ")")) %>%
  select(`Antiviral effectiveness scenario` = scenario_names, severity, estimate) %>%
  spread(severity, estimate) %>%
  arrange(`Antiviral effectiveness scenario`) 

tx_tab %>%
  kable(format = "latex", 
               align = "l|cc",
               caption = "\\textbf{Number of people undergoing prompt influenza antiviral treatment.} In these scenarios antiviral effectiveness against hospitalization is fixed at 20\\%.", 
               booktabs = TRUE) %>%
  row_spec(0:nrow(tx_tab), extra_latex_after = "\\addlinespace[3pt]")


